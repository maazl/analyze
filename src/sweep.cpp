#include "pca.h"
#include "parser.h"
#include "utils.h"
#include "pcmio.h"
#include "mathx.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <stdarg.h>

#include <chrono>
#include <thread>
#include <complex>
using namespace std;
typedef complex<double> Complex;

#define M_2PI (2.*M_PI)
#define M_3PI (3.*M_PI)
#ifndef M_PI_2
#define M_PI_2 (M_PI*.5)
#endif
#define M_PI_180 (M_PI/180.)
#define M_180_PI (180./M_PI)
#define M_sqrtPI (1.772453851)

#define SYNC_FIR 64


// config
static unsigned N = 65536;
static double freq = 48000;
static double f_min = 20;
static double f_max = 20000;
static double fstep = 1.05946309;
static unsigned discardsamp = 0;
static unsigned syncsamp = 20000;
static unsigned overlap = 1000;
static double synclevel = .3;
static double syncphase = 1;
static const char* infile = "-";
static const char* outfile = "-";
static const char* rawfile = NULL;
static const char* datafile = NULL;
static const char* zeroinfile = NULL; ///< file name for zero calibration data
static const char* zerooutfile = NULL;///< file name for differential zero calibration data
static const char* gaininfile = NULL; ///< file name for gain calibration data
static const char* gainoutfile = NULL;///< file name for differential gain calibration data
static unsigned loops = 1;
static bool verbose = false;
static bool floatin = false;
static bool swapbytes = false;
static const char* execcmd = NULL;

// data buffers
static unique_num_array<fftw_real> refbuffer;
static unique_num_array<fftw_real> inprebuffer1;
static unique_num_array<fftw_real> inbuffer1;
static unique_num_array<fftw_real> inprebuffer2;
static unique_num_array<fftw_real> inbuffer2;
static unique_num_array<double[2]> outprebuffer;
static unique_num_array<double[2]> outbuffer;

static FILE* resh = NULL;


// synchronize
static bool issync(double (*dp)[2])
{
	static int synccount = 0;
	double dphi = M_180_PI * fabs(fmod(dp[0][1] - dp[-SYNC_FIR][1] + M_2PI, M_2PI) - M_PI);
	if (verbose)
	{	static int cnt = 0;
		fprintf(stderr, "# Res: %2x %- 7f %- 7f %i %i %i %i %- 7f %- 7f\n",
			++cnt % SYNC_FIR, dp[0][0], M_180_PI * dp[0][1],
			dp[0][0] >= synclevel, dp[-SYNC_FIR][0] >= synclevel, dphi <= syncphase, synccount,
			dp[-SYNC_FIR][0], dphi);
	}

	if (dp[0][0] >= synclevel && dp[-SYNC_FIR][0] >= synclevel && dphi <= syncphase)
		return ++synccount == SYNC_FIR >> 2;
	synccount = 0;
	return false;
}

static int synchronize()
{  // FIR Filter
	static int rem = 0;
	const fftw_real* sp = inbuffer1.begin();
	const fftw_real* se = sp + (syncsamp >> 2);
	double (*dp)[2] = outbuffer.begin();
	sp += rem;
	while (sp < se)
	{
		const fftw_real* sp2 = sp;
		Complex x = 0;
		for (int l = SYNC_FIR >> 2; l; l--)
		{	x += Complex(sp2[0] - sp2[-2], sp2[-1] - sp2[-3]);
			sp2 -= 4;
		}
		(*dp)[0] = abs(x);
		(*dp)[1] = arg(x);
		//fprintf(stderr, "# Data: %4.4x %4.4x %4.4x %4.4x\t", fromraw(sp[0]), fromraw(sp[-2]), fromraw(sp[-4]), fromraw(sp[-6]));
		if (issync(dp))
			return (sp - inbuffer1.begin()) >> 1;
		dp += 2;
		sp += 8;
	}
	rem = sp - se;
	// save history
	memcpy(inprebuffer1.begin(), sp - SYNC_FIR, SYNC_FIR * sizeof(*inprebuffer1.begin()));
	memcpy(outprebuffer.begin(), dp - SYNC_FIR, SYNC_FIR * sizeof(*outprebuffer.begin()));
	return -1;
}

// analysis
struct ana
{	Complex Xi;
	double Sum;
	double Sum2;
	double noise() const { return sqrt(Sum2 - norm(Xi) - sqr(Sum)); }
	void store(Complex vec, double v)
	{	Xi += vec * v;
		Sum += v;
		Sum2 += sqr(v);
	}
	void finish()
	{	Xi /= N / M_SQRT2;
		Sum /= N;
		Sum2 /= N;
	}
} Ana[2];
Complex H;

static void analyze(int fi)
{
	memset(Ana, 0, sizeof Ana);
	const double fs = M_2PI * fi / N;
	const fftw_real* sp1 = inbuffer1.begin();
	const fftw_real* sp2 = inbuffer1.begin();
	for (unsigned i = 0; i < N; ++i)
	{	// calculate sin/cos sums
		auto phasevec = polar(1., fs * i);
		Ana[0].store(phasevec, *sp1++);
		Ana[1].store(phasevec, *sp2++);
	}
	Ana[0].finish();
	Ana[1].finish();
	H = Ana[1].Xi / Ana[0].Xi;
	//printf("X: %f\t%f\t%f\t%f\n", ana[6], ana[7], sqrt(ana[8]), sqrt(ana[9]));;
}

static void doanalysis()
{
	FILEguard in(infile, "rb");

	PCMinput pcmin(floatin ? Format::F32 : swapbytes ? Format::I16_SWAP : Format::I16);

	// skip initial samples
	pcmin.discard(in, discardsamp);

	// synchronize
	{	int i;
		unsigned n = 0;
		for (;;)
		{	if (n >= syncsamp)
				die(29, "Failed to syncronize.");
			pcmin.read(in, inbuffer1.slice(0, syncsamp>>2), inbuffer2.slice(0, syncsamp>>2), false);
			/*FILE* fs = fopen("sync.dat", "a");
			 write2ch(fs, inbuffer, syncsamp/4);
			 fclose(fs);*/
			i = synchronize();
			if (i >= 0)
				break;
			n += syncsamp / 4;
		}
		//fprintf(stderr, "Sync: %i\t%i\t%i\t%i\n", n, i, syncsamp/4, overlap);
		if ((int)(syncsamp >> 2) + i < SYNC_FIR)
			die(29, "Syncpoint missed by %i samples.", -((syncsamp >> 2) + i - SYNC_FIR));
		// synced. From now no samples must get lost.
		fprintf(stderr, "Synced after %i samples. Read %i samples, Skip another %i samples\n", n + i, n + (syncsamp >> 2), (syncsamp >> 2) + i - SYNC_FIR);
		// discard until end of sync - overlap
		pcmin.discard(in, (syncsamp >> 2) + i - SYNC_FIR);
	}
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
		{	// prepare write raw data
			FILEguard fr = NULL;
			if (rawfile)
			{	unique_ptr<char[]> buf(new char[strlen(rawfile + 20)]);
				sprintf(buf.get(), rawfile, findex * freq / N);
				fr = checkedopen(buf.get(), "w");
			}
			// discard first samples
			pcmin.read(in, inbuffer1, inbuffer2, false, fr);
			// show
			fputs("start...", stderr);
			pcmin.reset();
			// read data
			pcmin.read(in, inbuffer1, inbuffer2, false, fr);
			// show
			fputs("completed ", stderr);
		}
		// analyze data
		analyze(findex);
		// show
		fprintf(stderr, "%f.1 dB\n", 20 * log10(abs(H)));
		// write result
		fprintf(resh, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", findex * freq / N,
			abs(Ana[0].Xi), M_180_PI * arg(Ana[0].Xi), abs(Ana[1].Xi), M_180_PI * arg(Ana[1].Xi), abs(H), M_180_PI * arg(H),
			pcmin.Limits[0].Min, pcmin.Limits[0].Max, pcmin.Limits[1].Min, pcmin.Limits[0].Max, Ana[0].noise(), Ana[1].noise());
		fflush(resh);
		// next frequency
		fq *= fstep;
	}

}

// reference signal output
static void gensync()
{
	fftw_real* dp = refbuffer.begin();
	fftw_real* de = dp + syncsamp;
	while (dp != de)
	{
		dp[0] = -1;
		dp[1] = 0;
		dp[2] = 1;
		dp[3] = 0;
		dp += 4;
	}
	// phase jump: pi
	de += syncsamp;
	while (dp != de)
	{
		dp[0] = 1;
		dp[1] = 0;
		dp[2] = -1;
		dp[3] = 0;
		dp += 4;
	}
}

static void genref(int fi)
{  // fill reference buffer
	double fs = M_2PI * fi / N;
	for (unsigned i = 0; i < 2 * N; ++i)
		refbuffer[i] = cos(fs * i);
}

static void refplay()
{	// we should increase the priority here

	FILEguard out(outfile, "wb");

	PCMoutput pcmout(floatin ? Format::F32 : swapbytes ? Format::I16_SWAP : Format::I16);
	refbuffer.reset(N);

	// write wav header
	// TODO: calc exact length
	pcmout.WAVheader(out, ~0, freq);

	// pregap
	refbuffer.clear();
	pcmout.zero(out, discardsamp);
	// synchronize
	gensync();
	pcmout.write(out, refbuffer.slice(0, 2*syncsamp));
	// overlap
	pcmout.zero(out, overlap);
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
		pcmout.write(out, refbuffer);
		// write reference 2nd try
		pcmout.write(out, refbuffer);
		// next frequency
		fq *= fstep;
	}
}

/// Table of configuration parameters - MUST BE ORDERED BY NAME!
const reference<const OptionDesc> OptionMap[] =
{	MkOpt("bn",   "FFT lenght", N)
,	MkOpt("exec", "execute shell command", execcmd)
,	MkDOp("ff32", "32 bit floating point format", floatin, true)
,	MkDOp("fi16", "16 bit integer format (default)", floatin, false)
,	MkOpt("flog", "frequency increment factor", fstep)
,	MkOpt("fmax", "maximum frequency", f_max)
,	MkOpt("fmin", "minimum frequency", f_min)
,	MkOpt("fsamp","sampling rate", freq)
,	MkDOp("gg",   "generate gain calibration file", gainoutfile, "gain.dat")
,	MkDOp("gr",   "use gain calibration file", gaininfile, "gain.dat")
,	MkOpt("in",   "input file, stdin by default", infile)
,	MkOpt("ln",   "number of loops", loops)
,	MkSet("loop", "infinite mode", loops, UINT_MAX)
,	MkSet("ma",   "analyze only", outfile, (const char*)nullptr)
,	MkSet("mr",   "reference only", infile, (const char*)nullptr)
,	MkOpt("out",  "output file, stdout by default", outfile)
,	MkOpt("psa",  "discard first samples", discardsamp)
,	MkOpt("slvl", "sync level", synclevel)
,	MkOpt("sov",  "overlap", overlap)
,	MkOpt("sph",  "sync phase", syncphase)
,	MkOpt("sync", "sync samples", syncsamp)
,	MkOpt("v",    "verbose", verbose)
,	MkDOp("wd",   "write data file", datafile, "data.dat")
,	MkDOp("wr",   "write raw data file", rawfile, "raw.dat")
,	MkOpt("xb" ,  "swap bytes", swapbytes)
,	MkDOp("zg",   "generate matrix calibration file", zerooutfile, "zero.dat")
,	MkDOp("zr",   "use matrix calibration file", zeroinfile, "zero.dat")
};

int main(int argc, char* argv[])
{
	// parse cmdl
	{	Parser parser(OptionMap);
		while (--argc)
			parser.HandleArg(*++argv);
	}

	syncsamp &= ~7; // must be multiple of 8

	// setup input
	//PCMinput pcmin(floatin ? Format::F32 : swapbytes ? Format::I16_SWAP : Format::I16);
	// allocate buffers
	inprebuffer1.reset(N + SYNC_FIR);
	inprebuffer2.reset(N + SYNC_FIR);
	outprebuffer.reset(N + SYNC_FIR);
	inbuffer1.swap(inprebuffer1.slice(SYNC_FIR, N));
	inbuffer2.swap(inprebuffer2.slice(SYNC_FIR, N));
	outbuffer.swap(outprebuffer.slice(SYNC_FIR, N));

	inprebuffer1.clear();
	inprebuffer2.clear();
	outprebuffer.clear();

	resh = stdout;

	if (infile)
	{	if (outfile)
		{	// start reference generator
			thread(refplay).detach();
			this_thread::sleep_for(chrono::seconds(2));
		}
		// analyze
		doanalysis();
	} else if (outfile)
	{	// generate
		refplay();
	}

	// end reference player
	termrq = true;

	return 0;
}

