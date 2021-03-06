#include "utils.h"
#include "parser.h"
#include "mathx.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <ctype.h>


static double square(double p)
{
	return p < .5 ? 1. : -1.;
}

static double triangle(double p)
{
	return p < .5 ? 4. * p - 1. : 3. - 4. * p;
}

static double parabolic(double p)
{
	return p < .5 ? 1 - sqr(4. * p - 1.) : sqr(4. * p - 3.) - 1;
}

static double sine(double p)
{
	return sin(2 * M_PI * p);
}

static double saw_up(double p)
{
	return 2 * p - 1;
}

static double saw_down(double p)
{
	return 1 - 2 * p;
}

static double pulse(double p)
{
	return p == 0 ? 1 : -1;
}

static double (*parseshapefunc(const char* name))(double)
{
	static const struct fne
	{	char name[10];
		double (*fn)(double);
	} fnt[] =
	{	{ "parabolic", &parabolic },
		{ "pulse", &pulse },
		{ "saw-up", &saw_up },
		{ "saw-down", &saw_down },
		{ "sine", &sine },
		{ "square", &square },
		{ "triangle", &triangle } };

	if (name == NULL || *name == 0)
		return NULL;
	const fne* sp = (const fne*)bsearch(name, fnt, sizeof(fnt) / sizeof(*fnt), sizeof(*fnt), (int (*)(const void*, const void*))&strcasecmp);
	if (sp == NULL)
		die (34, "invalid shape %s", name);
	return sp->fn;
}

int main(int argc, char**argv)
{
	srand(clock());

	// command line
	const char* dfile = NULL;
	size_t nrep = 0;
	switch (argc)
	{default:
		die(45, "usage: %s nsamp[/harmonic] sampfreq shape [nrep [file]]]\n"
				"shape is one of square, triangle, sine, saw-up, saw-down, pulse, parabolic.\n", argv[0]);
	 case 6:
		dfile = argv[5];
	 case 5:
		nrep = atol(argv[4]);
	 case 4:;
	}
	size_t nsamp = atol(strtok(argv[1], "/"));
	const char* cp = strtok(NULL, "");
	size_t harmonic = 1;
	if (cp != NULL)
		harmonic = atol(cp);
	size_t sfreq = atol(argv[2]);
	double (*shapefn)(double) = parseshapefunc(argv[3]);

	// generate function
	short* buf = new short[2 * nsamp];
	short* dp = buf;
	for (size_t i = 0; i < nsamp; ++i)
	{
		short s = (short)floor(32767. * (*shapefn)(fmod((double)i / nsamp * harmonic, 1)) + myrand());
		*dp++ = s;
		*dp++ = -s;
	}

	// write
	FILEguard of = NULL;
	if (dfile)
	{	of = checkedopen(dfile, "wb");
		wavheader(of, 2 * nsamp * nrep, sfreq);
	} else // stdout
	{	of = binmode(stdout);
		wavheader(of, 0x1fffffdc, sfreq);
	}
	do
		fwrite(buf, 2 * nsamp * sizeof(short), 1, of);
	while (--nrep);

	// cleanup, well not in case of an exception...
	delete[] buf;
}

