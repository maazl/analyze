#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <stdarg.h>
#include <stdlib.h>

#include <complex>
using namespace std;
typedef complex<double> Complex;

#include "parser.h"

#define M_2PI (2.*M_PI)
#define M_3PI (3.*M_PI)
#ifndef M_PI_2
#define M_PI_2 (M_PI*.5)
#endif
#define M_PI_180 (M_PI/180.)
#define M_180_PI (180./M_PI)

#define N_MAX (65536*8)
#define CA_MAX 2

// data buffers
static short inbuffertmp[2 * CA_MAX * N_MAX];
static int sumbuffer[CA_MAX][65536];
static long nsamp;

void die(const char* msg, ...)
{
	va_list va;
	va_start(va, msg);
	vprintf(msg, va);
	va_end(va);
	putchar('\n');
	exit(1);
}

// config
static int N = 8192;
static bool swapbytes = false;
static bool writeraw = false;
static bool writedata = false;
static int discsamp = 0;
static int loops = 1;
static int binsz = 1;
static const char* infile = NULL;
static const char* execcmd = NULL;
static const char* plotcmd = NULL;

static inline short fromraw(short v)
{
	//return _srotl(v, 8);
	return (unsigned short)v >> 8 | v << 8;
}

static void asshort2(const short* src, size_t len)
{
	while (len--)
	{
		++sumbuffer[0][src[0] - SHRT_MIN];
		++sumbuffer[1][src[1] - SHRT_MIN];
		src += 2;
	}
}

static void asshortx2(const short* src, size_t len)
{
	while (len--)
	{
		++sumbuffer[0][fromraw(src[0]) - SHRT_MIN];
		++sumbuffer[1][fromraw(src[1]) - SHRT_MIN];
		src += 2;
	}
}

static inline double sqr(double v)
{
	return v * v;
}

static inline int64_t sqr(int64_t v)
{
	return v * v;
}

static inline double stddev(int64_t qsum, int64_t sum)
{  //printf("qs=%f s=%lli n=%li\n", (double)qsum, sum, nsamp);
	return sqrt((double)qsum / nsamp - sqr((double)sum / nsamp));
}

static void outdata()
{
	FILE* tout = fopen("data.dat", "wt");
	if (tout == NULL)
		die("Failed to create data.dat.");

	int* sp1 = sumbuffer[0];
	int* sp2 = sumbuffer[1];
	int64_t samp1 = 0;
	int64_t samp2 = 0;
	int64_t sum1 = 0;
	int64_t sum2 = 0;
	int64_t qsum1 = 0;
	int64_t qsum2 = 0;
	for (int s = SHRT_MIN; s <= SHRT_MAX; ++s, ++sp1, ++sp2)
	{  // write
		fprintf(tout, "%7i %12g %12g %7i %7i\n",
		// n  hl[n]                 hr[n]                 Nl[n] Nr[n]
		    s, (double)*sp1 / nsamp, (double)*sp2 / nsamp, *sp1, *sp2);
		samp1 += *sp1;
		samp2 += *sp2;
		int64_t v = (int64_t)*sp1 * s;
		sum1 += v;
		qsum1 += v * s;
		v = (int64_t)*sp2 * s;
		sum2 += v;
		qsum2 += v * s;
	}
	fclose(tout);
	if (samp1 != nsamp || samp2 != nsamp)
		die("Internal error: number of samples inconsistent. %li %li %li", samp1, samp2, nsamp);

	fprintf(stderr, "avg:\t%.2f\t%.2f\nstd:\t%.2f\t%.2f\n", (double)sum1 / nsamp, (double)sum2 / nsamp, stddev(qsum1, sum1), stddev(qsum2, sum2));
}

static void write1ch(FILE* out, const float* data, size_t len)
{
	while (len--)
		fprintf(out, "%g\n", *data++);
}

static void write2ch(FILE* out, const short* data, size_t len)
{
	while (len--)
	{
		fprintf(out, "%i\t%i\n", fromraw(data[0]), fromraw(data[1]));
		data += 2;
	}
}

static void write2ch(FILE* out, const float* data, size_t len)
{
	while (len--)
	{
		fprintf(out, "%g\t%g\n", data[0], data[1]);
		data += 2;
	}
}

static void write2ch(FILE* out, const float* data1, const float* data2, size_t len)
{
	while (len--)
		fprintf(out, "%g\t%g\n", *data1++, *data2++);
}

/*static char* abbrev(const char* s, const char* token)
 {  size_t l = strlen(s);
 return strnicmp(s, token, l) == 0 ? (char*)token + l : NULL;
 }*/

static void readN(const char* s, int* r)
{	bool ex;
	if ((ex = *s == '^'))
	++s;
	readint(s, r);
	if (ex)
	N = 1 << N;
}

const ArgMap argmap[] = // must be sorted
{	{ "bin",  (ArgFn)&readint, &binsz, 0}
,	{ "exec", (ArgFn)&readstring, &execcmd, 0}
,	{ "in" ,  (ArgFn)&readstring, &infile, 0}
,	{ "ln" ,  (ArgFn)&readint, &loops, 1}
,	{ "loop", (ArgFn)&setint, &loops, INT_MAX}
,	{ "n" ,   (ArgFn)&readN, &N, 0}
,	{ "plot", (ArgFn)&readstring, &plotcmd, 0}
,	{ "psa",  (ArgFn)&readintdef, &discsamp, 1}
,	{ "wd" ,  (ArgFn)&setflag, &writedata, true}
,	{ "wr" ,  (ArgFn)&setflag, &writeraw, true}
,	{ "xb" ,  (ArgFn)&setflag, &swapbytes, true}
};
const size_t argmap_size = sizeof argmap / sizeof *argmap;

int main(int argc, char* argv[])
{  // parse cmdl
	while (--argc)
		parsearg(*++argv);

	if (N > N_MAX)
		die("Data Length too large.");

	FILE* in;
	if (infile == NULL)
	{  // streaming
		in = stdin;
		//_fsetmode(in, "b");
	} else
	{
		in = fopen(infile, "rb");
		if (in == NULL)
			die("Failed to open input file.");
	}

	// discard first samples
	fread(inbuffertmp, sizeof *inbuffertmp, discsamp, in);

	memset(sumbuffer, 0, sizeof sumbuffer);
	nsamp = 0;

	// operation loop
	int loop = loops;
	do
	{
		if (fread(inbuffertmp, 2 * sizeof *inbuffertmp, N, in) != N)
			die("failed to read data");

		// write raw data
		if (writeraw)
		{
			FILE* tout = fopen("raw.dat", "wt");
			write2ch(tout, inbuffertmp, N);
			fclose(tout);
		}

		if (swapbytes)
			asshortx2(inbuffertmp, N);
		else
			asshort2(inbuffertmp, N);
		nsamp += N;

		// min/max analysis
		int minmax[2][2] = { { INT_MAX, INT_MIN }, { INT_MAX, INT_MIN } };
		int* sp1 = sumbuffer[0];
		for (int s = SHRT_MIN; s <= SHRT_MAX; ++s)
			if (*sp1++ != 0)
			{
				minmax[0][0] = s;
				break;
			}
		sp1 = sumbuffer[0] + 65536;
		for (int s = SHRT_MAX; s >= SHRT_MIN; --s)
			if (*--sp1 != 0)
			{
				minmax[0][1] = s;
				break;
			}
		int* sp2 = sumbuffer[1];
		for (int s = SHRT_MIN; s <= SHRT_MAX; ++s)
			if (*sp2++ != 0)
			{
				minmax[1][0] = s;
				break;
			}
		sp2 = sumbuffer[0] + 65536;
		for (int s = SHRT_MAX; s >= SHRT_MIN; --s)
			if (*--sp2 != 0)
			{
				minmax[1][1] = s;
				break;
			}
		// write raw status
		fprintf(stderr, "\nmin:\t%i\t%i\nmax:\t%i\t%i\n", minmax[0][0], minmax[1][0], minmax[0][1], minmax[1][1]);
		// write raw data
		/*if (writeraw)
		 {  tout = fopen("raw.dat", "wt");
		 write2ch(tout, inbuffer1, inbuffer2, N);
		 fclose(tout);
		 }*/

		if (writedata)
			outdata();

		if (plotcmd)
		{  // for gnuplot!
			puts(plotcmd);
			fflush(stdout);
		}

	} while (--loop);
	// close stdin to signal playrec
	fclose(in);

	outdata();

	return 0;
}
