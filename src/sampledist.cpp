#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>

#include "parser.h"
#include "utils.h"
#include "moment.h"

#define M_2PI (2.*M_PI)
#define M_3PI (3.*M_PI)
#ifndef M_PI_2
#define M_PI_2 (M_PI*.5)
#endif
#define M_PI_180 (M_PI/180.)
#define M_180_PI (180./M_PI)

#define N_MAX (65536*8)

// data buffers
static scoped_array<short> inbuffertmp;
static int sumbuffer[2][65536];
static int64_t nsamp;

// config
static unsigned N = 32768;
static bool swapbytes = false;
static const char* rawfile = NULL; // file name for raw data
static const char* datafile = NULL; // file name for histogram data
static unsigned discsamp = 0;
static unsigned addloop = UINT_MAX;
static unsigned loops = 1;
static const char* infile = NULL;
static const char* execcmd = NULL;
static const char* plotcmd = NULL;

static void asshort2(const short* src, size_t len)
{
	while (len--)
	{	++sumbuffer[0][src[0] - SHRT_MIN];
		++sumbuffer[1][src[1] - SHRT_MIN];
		src += 2;
	}
}

static void asshortx2(const short* src, size_t len)
{
	while (len--)
	{	++sumbuffer[0][bswap(src[0]) - SHRT_MIN];
		++sumbuffer[1][bswap(src[1]) - SHRT_MIN];
		src += 2;
	}
}

static void analyze()
{
	// min/max analysis
	int minmax[2][2] = { { INT_MAX, INT_MIN }, { INT_MAX, INT_MIN } };
	int* sp = sumbuffer[0];
	for (int s = SHRT_MIN; s <= SHRT_MAX; ++s)
		if (*sp++ != 0)
		{	minmax[0][0] = s;
			break;
		}
	sp = sumbuffer[0] + 65536;
	for (int s = SHRT_MAX; s >= SHRT_MIN; --s)
		if (*--sp != 0)
		{	minmax[0][1] = s;
			break;
		}
	sp = sumbuffer[1];
	for (int s = SHRT_MIN; s <= SHRT_MAX; ++s)
		if (*sp++ != 0)
		{	minmax[1][0] = s;
			break;
		}
	sp = sumbuffer[1] + 65536;
	for (int s = SHRT_MAX; s >= SHRT_MIN; --s)
		if (*--sp != 0)
		{	minmax[1][1] = s;
			break;
		}
	// write raw status
	fputs("\n\tsamples \tdB\n", stderr);
	fprintf(stderr, "min:\t%i\t%i\t%.1f\t%.1f\n", minmax[0][0], minmax[1][0], todB(-minmax[0][0] / 32767.), todB(-minmax[1][0] / 32767.));
	fprintf(stderr, "max:\t%i\t%i\t%.1f\t%.1f\n", minmax[0][1], minmax[1][1], todB(minmax[0][1] / 32767.), todB(minmax[1][1] / 32767.));

	int* sp1 = sumbuffer[0];
	int* sp2 = sumbuffer[1];
	Kurtosis stat1;
	Kurtosis stat2;
	for (int s = SHRT_MIN; s <= SHRT_MAX; ++s)
	{	// write
		stat1.push(s, *sp1++);
		stat2.push(s, *sp2++);
	}
	assert(stat1.count() == nsamp && stat2.count() == nsamp);

	fprintf(stderr, "mean\t%.2f\t%.2f\n", stat1.mean(), stat2.mean());
	fprintf(stderr, "stddev\t%.1f\t%.1f\n", stat1.stddeviation(), stat2.stddeviation());
	fprintf(stderr, "skew\t%.4f\t%.4f\n", stat1.skewness(), stat2.skewnessBC());
	fprintf(stderr, "kurtos.\t%.4f\t%.4f\n", stat1.kurtosisExcessBC(), stat2.kurtosisExcessBC());
	double crest1 = (minmax[0][1] - minmax[0][0]) / sqrt(stat1.sum2() / stat1.count()) / 2;
	double crest2 = (minmax[1][1] - minmax[1][0]) / sqrt(stat2.sum2() / stat2.count()) / 2;
	fprintf(stderr, "crest\t%.4f\t%.4f\t%.1f\t%.1f\n", crest1, crest2, todB(crest1), todB(crest2));
}

static void outdata()
{	FILEguard tout(datafile, "wt");
	int* sp1 = sumbuffer[0];
	int* sp2 = sumbuffer[1];
	fputs("#n\thl\thr\tNl\tNr\n", tout);
	for (int s = SHRT_MIN; s <= SHRT_MAX; ++s, sp1++, sp2++)
		fprintf(tout, "%7i\t%12g\t%12g\t%7i\t%7i\n",
		// n  hl[n]                 hr[n]                 Nl[n] Nr[n]
		   s, (double)*sp1 / nsamp, (double)*sp2 / nsamp, *sp1, *sp2);
}

static void write2ch(FILE* out, const short* data, size_t len)
{
	while (len--)
	{	fprintf(out, "%i\t%i\n", bswap(data[0]), bswap(data[1]));
		data += 2;
	}
}

const OptionDesc OptionMap[] = // must be sorted
{	MkOpt("al",   "add blocks, infinite by default", &addloop)
,	MkOpt("bn",   "block length, 32768 by default", &N)
,	MkOpt("df",   "write histogram data to file", &datafile)
,	MkOpt("exec", "execute shell command after block completed", &execcmd)
,	MkOpt("in",   "name of input file, default stdin", &infile)
,	MkOpt("ln",   "number of cycles, 1 by default", &loops)
,	MkOpt("loop", "infinite input", &loops, UINT_MAX)
,	MkOpt("plot", "pipe command after block completed", &plotcmd)
,	MkOpt("psa",  "discard first samples", &discsamp)
,	MkOpt("rf",   "write raw data file (diagnostics)", &rawfile)
,	MkOpt("wd",   "write histogram data to hist.dat", &datafile, "hist.dat")
,	MkOpt("wr",   "write raw data to raw.dat", &rawfile, "raw.dat")
,	MkOpt("xb" ,  "swap bytes", &swapbytes)
};

int main(int argc, char* argv[])
{
	// parse cmdl
	{	Parser parser(OptionMap);
		while (--argc)
			parser.HandleArg(*++argv);
	}

	if (N > N_MAX)
		die(32, "Data Length too large.");

	// allocate buffers
	inbuffertmp.reset(2 * N);

	FILEguard in = infile == NULL
		? binmode(stdin)
		: checkedopen(infile, "rb");

	// discard first samples
	fread(inbuffertmp.begin(), sizeof(short), 2 * discsamp, in);

	memset(sumbuffer, 0, sizeof sumbuffer);
	nsamp = 0;

	// operation loop
	unsigned addloops = 0;
	unsigned loop = loops;
	do
	{
		fread2(inbuffertmp.begin(), 2 * sizeof(short), N, in);

		// write raw data
		if (rawfile)
			write2ch(FILEguard(rawfile, "wt"), inbuffertmp.begin(), N);

		if (swapbytes)
			asshortx2(inbuffertmp.begin(), N);
		else
			asshort2(inbuffertmp.begin(), N);
		nsamp += N;

		analyze();

		if (datafile)
			outdata();

		if (execcmd)
			system(execcmd);
		if (plotcmd)
		{	// for gnuplot!
			puts(plotcmd);
			fflush(stdout);
		}

		if (++addloops < addloop)
			continue; // add more data
		memset(sumbuffer, 0, sizeof sumbuffer);
		nsamp = 0;
		addloops = 0;

	} while (--loop);

	return 0;
}
