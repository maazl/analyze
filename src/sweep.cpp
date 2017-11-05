#include "pca.h"
#include "parser.h"
#include "utils.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <stdarg.h>

#ifdef __OS2__
#define INCL_BASE
#include <os2.h>
#endif

#include <chrono>
#include <thread>
#include <complex>
#include <fftw3.h>
using namespace std;
typedef double fftw_real;
typedef complex<double> Complex;

#define M_2PI (2.*M_PI)
#define M_3PI (3.*M_PI)
#ifndef M_PI_2
#define M_PI_2 (M_PI*.5)
#endif
#define M_PI_180 (M_PI/180.)
#define M_180_PI (180./M_PI)
#define M_sqrtPI (1.772453851)

#define N_MAX (65536*8)
#define CA_MAX 2
#define SYNC_FIR 64

int __gxx_personality_v0; // gcc @�$%&!

// data buffers
static short refbuffer[N_MAX * 2];
static short inprebuffer[CA_MAX * (N_MAX + SYNC_FIR)];
static short* const inbuffer = inprebuffer + CA_MAX * SYNC_FIR;
static float outprebuffer[N_MAX + 2 * SYNC_FIR];
static float* const outbuffer = outprebuffer + 2 * SYNC_FIR;

static FILE* resh = NULL;

static const double minval = 1E-20;

static void vectorscale(fftw_real* data, double factor, size_t len)
{
	while (len--)
		*data++ *= factor;
}

static void vectoradd(fftw_real* dst, fftw_real* src, size_t len)
{
	while (len--)
		*dst++ += *src++;
}

static void vektormul(fftw_real* dst, fftw_real* src, size_t len)
{
	while (len--)
		*dst++ *= *src++;
}

// config
static unsigned N = 8192;
static double freq = 48000;
static double f_min = 20;
static double f_max = 20000;
static double fstep = 1.05946309;
static unsigned discardsamp = 0;
static unsigned syncsamp = 50000;
static unsigned overlap = 1000;
static double synclevel = 10000;
static double syncphase = 1;
static bool writeraw = false;
static bool writedata = false;
static unsigned loops = 1;
static int mode = 3; // 1 = generate, 2 = analyze 3 = both
static int zeromode = 0; // 0 = none, 1 = read, 2 = generate
static int gainmode = 0; // 0 = none, 1 = read, 2 = generate, 3 = generatedelta
static bool verbose = false;
static const char* infile = NULL;
static const char* execcmd = NULL;

static inline short fromraw(short v)
{
	//return _srotl(v, 8);
	return (unsigned short)v >> 8 | v << 8;
}

static int minmax[4];

static inline short storeminmax(short val, int* dst)
{
	if (val < dst[0])
		dst[0] = val;
	if (val > dst[1])
		dst[1] = val;
	return val;
}

static void short2float2(fftw_real* dst1, fftw_real* dst2, const short* src, size_t len)
{
	while (len)
	{
		*dst1++ = storeminmax(fromraw(src[0]), minmax);
		*dst2++ = storeminmax(fromraw(src[1]), minmax + 2);
		--len;
		src += 2;
	}
}

static inline double abs(double d1, double d2)
{
	return sqrt(sqr(d1) + sqr(d2));
}

void init()
{
	minmax[0] = INT_MAX;
	minmax[1] = INT_MIN;
	minmax[2] = INT_MAX;
	minmax[3] = INT_MIN;
}

/*void write1ch(FILE* out, const float* data, size_t len)
 {  while (len--)
 fprintf(out, "%g\n", *data++);
 }*/

void write2ch(FILE* out, const short* data, size_t len)
{
	while (len--)
	{
		fprintf(out, "%i\t%i\n", fromraw(data[0]), fromraw(data[1]));
		data += 2;
	}
}

/*void write2ch(FILE* out, const float* data, size_t len)
 {  while (len--)
 {  fprintf(out, "%g\t%g\n", data[0], data[1]);
 data += 2;
 }
 }

 void write2ch(FILE* out, const float* data1, const float* data2, size_t len)
 {  while (len--)
 fprintf(out, "%g\t%g\n", *data1++, *data2++);
 }

 void writepolar(FILE* out, const float* data, size_t len, double inc)
 {  const float* data2 = data + len;
 // 1st line
 fprintf(out, "0\t%g\t0\n", *data++);
 len = 1;
 while (data < --data2)
 fprintf(out, "%g\t%g\t%g\n", len++*inc, *data++, *data2);
 }

 void writecomplex(FILE* out, const Complex* data, size_t len)
 {  while (len--)
 {  fprintf(out, "%14g\t%14g\t%14g\t%14g\n", data->real(), data->imag(), abs(*data), arg(*data)*M_180_PI);
 ++data;
 }  }

 void readcomplex(FILE*in, Complex* data, size_t len)
 {  while (len--)
 {  double a,b;
 if (fscanf(in, "%lg%lg%*[^\n]", &a, &b) != 2)
 die("Failed to read complex data (%i).", errno);
 //(stderr, "%g\t%g\n", a,b);
 *data++ = Complex(a,b);
 }  }*/

static void fwriteexact(const void* buffer, size_t size, size_t n, FILE* f)
{
	while (n)
	{
		size_t r = fwrite(buffer, size, n, f);
		if (r <= 0)
			die(27, "Failed to write %lu blocks a %lu bytes", n, size);
		n -= r;
		(const char*&)buffer += r;
	}
}

static void freadexact(void* buffer, size_t size, size_t n, FILE* f)
{
	while (n)
	{
		size_t r = fread(buffer, size, n, f);
		if (r <= 0)
			die(27, "Failed to read %lu blocks a %lu bytes", n, size);
		//fwrite(buffer, size, r, stdout);
		n -= r;
		(char*&)buffer += r;
	}
}

// synchronize
static void syncinit()
{
	memset(inprebuffer, 0, sizeof inprebuffer);
	memset(outprebuffer, 0, sizeof outprebuffer);
}

static bool issync(float* dp)
{
	static int synccount = 0;
	if (verbose)
	{
		static int cnt = 0;
		fprintf(stderr, "# Res: %2x %- 7f %- 7f %i %i %i %i %- 7f %- 7f %- 7f\n", ++cnt % SYNC_FIR, dp[0], M_180_PI * dp[1], dp[0] >= synclevel,
		    dp[-2 * SYNC_FIR] >= synclevel, M_180_PI * abs(fmod(dp[1] - dp[1 - 2 * SYNC_FIR] + M_2PI, M_2PI) - M_PI) <= syncphase, synccount, dp[-2 * SYNC_FIR],
		    M_180_PI * abs(fmod(dp[1] - dp[1 - 2 * SYNC_FIR] + M_2PI, M_2PI) - M_PI), M_180_PI * abs(fmod(dp[1] - dp[1 - 2 * SYNC_FIR] + M_2PI, M_2PI) - M_PI));
	}

	if (dp[0] >= synclevel && dp[-2 * SYNC_FIR] >= synclevel && M_180_PI * abs(fmod(dp[1] - dp[1 - 2 * SYNC_FIR] + M_2PI, M_2PI) - M_PI) <= syncphase)
		return ++synccount == SYNC_FIR >> 2;
	synccount = 0;
	return false;
}

static int synchronize(int len)
{  // FIR Filter
	static int rem = 0;
	const short* sp = inbuffer;
	float* dp = outbuffer;
	const short* se = sp + 2 * len;
	sp += rem;
	while (sp < se)
	{
		const short* sp2 = sp;
		int a = 0;
		int b = 0;
		for (int l = SYNC_FIR >> 2; l; l--)
		{
			a += fromraw(sp2[0]) - fromraw(sp2[-4]);
			b += fromraw(sp2[-2]) - fromraw(sp2[-6]);
			sp2 -= 8;
		}
		dp[0] = abs(a, b);
		dp[1] = atan2((float)b, a);
		//fprintf(stderr, "# Data: %4.4x %4.4x %4.4x %4.4x\t", fromraw(sp[0]), fromraw(sp[-2]), fromraw(sp[-4]), fromraw(sp[-6]));
		if (issync(dp))
			return (sp - inbuffer) >> 1;
		dp += 2;
		sp += 8;
	}
	rem = sp - se;
	// save history
	memcpy(inprebuffer, sp - 2 * SYNC_FIR, 2 * SYNC_FIR * sizeof(short));
	memcpy(outprebuffer, dp - 2 * SYNC_FIR, 2 * SYNC_FIR * sizeof(float));
	return -1;
}

// analysis
static double ana[10];

static void analyze(int fi)
{
	init();
	memset(ana, 0, sizeof ana);
	const double fs = M_2PI * fi / N;
	const short* sp = inbuffer;
	for (int i = 0; i < N; ++i)
	{
		double s = sin(fs * i);
		double c = cos(fs * i);
		// calculate sin/cos sums
		double v = storeminmax(fromraw(sp[0]), minmax);
		ana[0] += c * v;
		ana[1] += s * v;
		ana[6] += v;
		ana[8] += sqr(v);
		v = storeminmax(fromraw(sp[1]), minmax + 2);
		ana[2] += c * v;
		ana[3] += s * v;
		ana[7] += v;
		ana[9] += sqr(v);
		sp += 2;
	}
	ana[0] /= N / M_SQRT2;
	ana[1] /= N / M_SQRT2;
	ana[2] /= N / M_SQRT2;
	ana[3] /= N / M_SQRT2;
	ana[6] /= N;
	ana[7] /= N;
	ana[8] /= N;
	ana[9] /= N;
	double s = sqr(ana[0]) + sqr(ana[1]);
	ana[4] = (ana[0] * ana[2] + ana[1] * ana[3]) / s;
	ana[5] = (ana[0] * ana[3] - ana[1] * ana[2]) / s;
	//printf("X: %f\t%f\t%f\t%f\n", ana[6], ana[7], sqrt(ana[8]), sqrt(ana[9]));
	ana[6] = sqrt(ana[8] - s - sqr(ana[6]));
	ana[7] = sqrt(ana[9] - sqr(ana[2]) - sqr(ana[3]) - sqr(ana[7]));
}

static void doanalysis()
{
	FILEguard in = infile == NULL
		? binmode(stdin)
		: checkedopen(infile, "rb");

	// skip initial samples
	freadexact(inbuffer, 2 * sizeof(short), discardsamp, in);

	// synchronize
	int i;
	int n = 0;
	for (;;)
	{
		if (n >= syncsamp)
			die(29, "Failed to syncronize.");
		freadexact(inbuffer, 2 * sizeof(short), syncsamp / 4, in);
		/*FILE* fs = fopen("sync.dat", "a");
		 write2ch(fs, inbuffer, syncsamp/4);
		 fclose(fs);*/
		i = synchronize(syncsamp >> 2);
		if (i >= 0)
			break;
		n += syncsamp / 4;
	}
	//fprintf(stderr, "Sync: %i\t%i\t%i\t%i\n", n, i, syncsamp/4, overlap);
	if ((syncsamp >> 2) + i - SYNC_FIR < 0)
		die(29, "Syncpoint missed by %i samples.", -((syncsamp >> 2) + i - SYNC_FIR));
	// synced. From now no samples must get lost.
	fprintf(stderr, "Synced after %i samples. Read %i samples, Skip another %i samples\n", n + i, n + (syncsamp >> 2), (syncsamp >> 2) + i - SYNC_FIR);
	// discard until end of sync - overlap
	freadexact(inbuffer, 2 * sizeof(short), (syncsamp >> 2) + i - SYNC_FIR, in);

	// Scan !
	double fq = f_min;
	int findex = 0;
	while (fq < f_max)
	{
		++findex;
		int nfi = (int)(fq / freq * N + .5);
		if (findex < nfi)
			findex = nfi;
		// show
		fprintf(stderr, "Now at %.2f Hz ", findex * freq / N);
		// discard first samples
		freadexact(inbuffer, 2 * sizeof(short), N, in);
		// write raw data #1
		FILEguard fr = NULL;
		if (writeraw)
		{	char buf[30];
			sprintf(buf, "raw_%3f", findex * freq / N);
			buf[8] = 0;
			strcat(buf, ".dat");
			fr = checkedopen(buf, "w");
			write2ch(fr, inbuffer, N);
		}
		// show
		fputs("start...", stderr);
		// read data
		freadexact(inbuffer, 2 * sizeof(short), N, in);
		// write raw data
		if (writeraw)
		{
			write2ch(fr, inbuffer, N);
			fr = NULL;
		}
		// show
		fputs("completed ", stderr);
		// analyze data
		analyze(findex);
		// show
		fprintf(stderr, "%f.1 dB\n", 20 * log10(sqrt(sqr(ana[4]) + sqr(ana[5]))));
		// write result
		fprintf(resh, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%i\t%f\t%f\n", findex * freq / N, sqrt(sqr(ana[0]) + sqr(ana[1])), M_180_PI * atan2(ana[1], ana[0]),
		    sqrt(sqr(ana[2]) + sqr(ana[3])), M_180_PI * atan2(ana[3], ana[2]), sqrt(sqr(ana[4]) + sqr(ana[5])), M_180_PI * atan2(ana[5], ana[4]), minmax[0],
		    minmax[1], minmax[2], minmax[3], ana[6], ana[7]);
		fflush(resh);
		// next frequency
		fq *= fstep;
	}

}

// reference signal output
static void gensync()
{
	short* dp = refbuffer;
	short* de = dp + syncsamp;
	while (dp != de)
	{
		dp[0] = -32767;
		dp[1] = 32767;
		dp[2] = 0;
		dp[3] = 0;
		dp[4] = 32767;
		dp[5] = -32767;
		dp[6] = 0;
		dp[7] = 0;
		dp += 8;
	}
	// phase jump: Pi
	de = refbuffer + 2 * syncsamp;
	while (dp != de)
	{
		dp[0] = 32767;
		dp[1] = -32767;
		dp[2] = 0;
		dp[3] = 0;
		dp[4] = -32767;
		dp[5] = 32767;
		dp[6] = 0;
		dp[7] = 0;
		dp += 8;
	}
}

static void genref(int fi)
{  // fill reference buffer
	double fs = M_2PI * fi / N;
	for (int i = 0; i < 2 * N; i += 2)
		refbuffer[i + 1] = -(refbuffer[i] = (short)(32767 * cos(fs * i / 2) + rand() / (RAND_MAX + 1.)));
}

static void refplay()
{  // we should increase the priority here
#ifdef __OS2__
	DosSetPriority(PRTYS_THREAD, PRTYC_NOCHANGE, 1, 0);
#else
//#error You need to implement a platform dependent way to increase the priority of this thread.
#endif

	// write wav header
	//_fsetmode(stdout, "b");
	static const char wavhdr[44] =
	{	'R', 'I', 'F', 'F', '\x94', '\xff', '\xff', 0x7f,
		'W', 'A', 'V', 'E', 'f', 'm', 't', ' ',
		16, 0, 0, 0, 1, 0, 2, 0,
		'\x80', '\xbb', 0, 0, 0, '\xee', 2, 0, 4, 0, 16, 0,
		'd', 'a', 't', 'a', 0x70, '\xff', '\xff', 0x7f
	};
	fwriteexact(wavhdr, 1, sizeof wavhdr, stdout);

	// pregap
	memset(refbuffer, 0, discardsamp * 2 * sizeof(short));
	fwriteexact(refbuffer, 2 * sizeof(short), discardsamp, stdout);
	// synchronize
	gensync();
	fwriteexact(refbuffer, 2 * sizeof(short), syncsamp, stdout);
	// overlap
	memset(refbuffer, 0, overlap * 2 * sizeof(short));
	fwriteexact(refbuffer, 2 * sizeof(short), overlap, stdout);
	// Wobble !
	double fq = f_min;
	int findex = 0;
	while (fq < f_max && !termrq)
	{
		++findex;
		int nfi = (int)(fq / freq * N + .5);
		if (findex < nfi)
			findex = nfi;
		genref(findex);
		// write reference
		fwriteexact(refbuffer, 2 * sizeof(short), N, stdout);
		// write reference 2nd try
		fwriteexact(refbuffer, 2 * sizeof(short), N, stdout);
		// next frequency
		fq *= fstep;
	}
}

const OptionDesc OptionMap[] = // must be sorted
{	MkOpt("bn",   "FFT lenght", &N)
,	MkOpt("exec", "execute shell command", &execcmd)
,	MkOpt("flog", "frequency increment factor", &fstep)
,	MkOpt("fmax", "maximum frequency", &f_max)
,	MkOpt("fmin", "minimum frequency", &f_min)
,	MkOpt("fsamp","sampling rate", &freq)
,	MkOpt("gd",   "verify gain calibration", &gainmode, 3)
,	MkOpt("gg",   "generate gain calibration file", &gainmode, 2)
,	MkOpt("gr",   "use gain calibration file", &gainmode, 1)
,	MkOpt("in",   "input file, stdin by default", &infile)
,	MkOpt("ln",   "number of loops", &loops)
,	MkOpt("loop", "infinite mode", &loops, UINT_MAX)
,	MkOpt("ma",   "analyze only", &mode, 2)
,	MkOpt("mr",   "reference only", &mode, 1)
,	MkOpt("psa",  "discard first samples", &discardsamp)
,	MkOpt("slvl", "sync level", &synclevel)
,	MkOpt("sov",  "overlap", &overlap)
,	MkOpt("sph",  "sync phase", &syncphase)
,	MkOpt("sync", "sync samples", &syncsamp)
,	MkOpt("v",    "verbose", &verbose)
,	MkOpt("wd",   "write data", &writedata)
,	MkOpt("wr",   "write raw data", &writeraw)
,	MkOpt("zg",   "generate matrix calibration file", &zeromode, 2)
,	MkOpt("zr",   "use matrix calibration file", &zeromode, 1)
};

int main(int argc, char* argv[])
{
	// parse cmdl
	{	Parser parser(OptionMap);
		while (--argc)
			parser.HandleArg(*++argv);
	}

	syncsamp &= ~7; // must be multiple of 8

	/*resh = fdopen(3, "w");
	 if (resh == NULL)
	 die("Failed to open output handle 3 (%i).", errno);*/
	resh = stdout;

	switch (mode)
	{
	case 1:
		// generate
		refplay();
		break;
	case 3:
		{	// start reference generator
			thread(refplay).detach();
			this_thread::sleep_for(chrono::seconds(2));
		}
	case 2:
		// analyze
		doanalysis();
		break;
	default:
		die(32, "Unsupported mode: %i", mode);
	}

	// end reference player
	termrq = true;

	return 0;
}

