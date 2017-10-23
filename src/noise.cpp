#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#include <fftw3.h>

#include <complex>
using namespace std;
typedef double fftw_real;
typedef complex<double> Complex;

#include "parser.h"

#define M_2PI (2.*M_PI)
#define M_180_PI (180./M_PI)


static unsigned    n_fft      = 0;     // period in samples
static unsigned    f_samp     = 48000; // sampling rate
static double      mgain      = 0;     // master gain [dB]
static double      f_min      = 0;
static double      f_max      = 0;
static double      f_inc      = 1;
static double      f_log      = 0;
static double      scalepow   = 0;
static unsigned    n_harmonic = 0;
static bool        stereo     = false;
static bool        sweep      = false;
static size_t      synclen    = 0;
static const char* F_data     = NULL;
static const char* F_res      = NULL;
static const char* resmode    = "w";
static const char* F_wav      = NULL;
static unsigned    n_rep      = 1;
//static int chirpphase = 0;
static const char* execcmd    = NULL;  // shell command to execute after analysis

static double      mfact; // Master gain factor

void wavwriter(FILE* fo, size_t n_samp)
{	// fake wav header
	unsigned wavhdr[11] =
	// 48000
	{	0x46464952, 0xffffffff, 0x45564157, 0x20746D66,
		0x00000010, 0x00020001, 0x0000BB80, 0x0002ee00,
		0x00100004, 0x61746164, 0xffffffff };

	// patch data size
	wavhdr[10] = n_samp * sizeof(short);
	wavhdr[1] = wavhdr[10] + 44;
	// patch sampling rate
	wavhdr[6] = f_samp;
	wavhdr[7] = f_samp << 2; // * 4 bytes per sample

	fwrite(wavhdr, sizeof wavhdr, 1, fo);
}

inline static double myrand()
{	unsigned long l = (((rand() << 11) ^ rand()) << 11) ^ rand();
	return (double)l / ULONG_MAX;
}


static void readN(const char* s, unsigned* r)
{	bool ex;
	if ((ex = *s == '^'))
		++s;
	readuint(s, r);
	if (ex)
		*r = 1 << *r;
}

static void quantize1(short* dst, const fftw_real* src, size_t count)
{	const fftw_real* const spe = src + count;
	while (src != spe)
	{	register short s = (short)floor(*src++ * mfact + myrand());
		dst[0] = s;
		dst[1] = -s;
		dst += 2;
	}
}

static void quantize2(short* dst, const fftw_real* src1, const fftw_real* src2, size_t count)
{	const fftw_real* const spe = src1 + count;
	while (src1 != spe)
	{	dst[0] = (short)floor(*src1++ * mfact + myrand());
		dst[1] = (short)floor(*src2++ * mfact + myrand());
		dst += 2;
	}
}

const ArgMap argmap[] = // must be sorted
{	{"ar",   (ArgFn)&setint,     &resmode,     (long)"a"}
,	{"bn"  , (ArgFn)&readN,      &n_fft,       0}
,	{"exec", (ArgFn)&readstring, &execcmd,     0}
,	{"finc", (ArgFn)&readdouble, &f_inc,       0}
,	{"flog", (ArgFn)&readdouble, &f_log,       0}
,	{"fmax", (ArgFn)&readdouble, &f_max,       0}
,	{"fmin", (ArgFn)&readdouble, &f_min,       0}
,	{"fsamp",(ArgFn)&readuint,   &f_samp,      0}
,	{"gm",   (ArgFn)&readdouble, &mgain,       0}
,	{"harm", (ArgFn)&readuintdef,&n_harmonic,  3}
,	{"ln" ,  (ArgFn)&readuint,   &n_rep,       1}
,	{"loop", (ArgFn)&setuint,    &n_rep,       0}
,	{"mst",  (ArgFn)&setflag,    &stereo,      true}
,	{"msweep",(ArgFn)&setflag,   &sweep,       true}
,	{"scale",(ArgFn)&readdouble, &scalepow,    0}
,	{"sync", (ArgFn)&readintdef, &synclen,     2}
,	{"wd",   (ArgFn)&readstring, &F_data,      0}
,	{"wr",   (ArgFn)&readstring, &F_res,       0}
,	{"ww",   (ArgFn)&readstring, &F_wav,       0}
};
const size_t argmap_size = sizeof argmap / sizeof *argmap;


int main(int argc, char**argv)
{	srand(clock());

	// parse cmdl
	while (--argc)
		parsearg(*++argv);

	if (n_fft == 0)
		die(48, "usage: noise bn<fft_size> fmin<min_freq> fmax<max_freq> [wd<data_file>] [wr<wave_file>]\n"
				  "See documentation for more options.");

	mfact = 32767. * pow(10., mgain/20.);

	// round fmin & fmax
	double f_bin = (double)f_samp/n_fft;
	size_t i_min = (int)floor(f_min/f_bin +.5);
	size_t i_max = (int)floor(f_max/f_bin +.5);
	if (i_max > n_fft/2 || i_min > i_max)
		die(34, "fmin and/or fmax out of range");
	f_inc -= .5;
	f_log += 1;
	fprintf(stderr, "imin=%zi imax=%zi finc=%f flog=%f\n", i_min, i_max, f_inc, f_log);

	// generate coefficients
	fftw_real* fftbuf = fftw_alloc_real((stereo+1)*(n_fft+1));
	int* harmonics = new int[n_fft+2];
	if (fftbuf == NULL || harmonics == NULL)
		die(39, "malloc(%lu) failed", n_fft);
	memset(fftbuf, 0, (stereo+1)*(n_fft+1) * sizeof *fftbuf);   // all coefficients -> 0
	memset(harmonics, 0, (n_fft/2+1) * sizeof *harmonics);
	size_t fcount = 0; // number of used frequencies
	float maxamp = 0;
	int sign = 1;
	for (size_t i = i_min; i <= i_max; ++i)
	{	if (i)
		{	// skip used harmonics
			for (size_t j = 1; j <= n_harmonic && i*j <= n_fft/2; ++j)
				if (harmonics[i*j])
					goto next_f;
			// lock harmonics
			for (size_t j = 1; i*j <= n_fft/2; ++j)
				harmonics[i*j] = j * sign;
		}
		// calculate coefficients
		++fcount;
		fftbuf[i] = pow(i, scalepow);
		if (fftbuf[i] > maxamp)
			maxamp = fftbuf[i];
		if (!sweep && i && i != n_fft/2)
		{	// apply random phase
			double phi = 2*M_PI * myrand();
			fftbuf[n_fft-i] = fftbuf[i] * sin(phi); // b[i]
			fftbuf[i] *= cos(phi);                  // a[i]
		}
		// second channel
		if (sign < 0)
		{	fftbuf[2*n_fft+1-i] = fftbuf[n_fft-i];
			fftbuf[n_fft+1+i] = fftbuf[i];
		}
		// next frequency
		if (stereo && !sweep)
			sign = -sign;
		//fprintf(stderr, "f %i %i\n", i, (int)floor(i * f_log + f_inc));
		i = (int)floor(i * f_log + f_inc);
	 next_f:;
	}
	// normalize
	{	fftw_real* dp = fftbuf;
		fftw_real* ep = dp + (stereo+1)*(n_fft+1);
		while (dp != ep)
			*dp++ /= maxamp;
	}

	// write design result
	if (F_data)
	{	FILE* of = fopen(F_data, "w");
		if (of == NULL)
			die(41, "Failed to open %s for writing", F_data);
		for (size_t i = 0; i <= n_fft/2; ++i)
		{	Complex ci(fftbuf[i], i && i != n_fft/2 ? fftbuf[n_fft-i] : 0);
			//          freq  |A|  arg A ai  bi   harm
			fprintf(of, "%12g %12g %12g %12g %12g %8i\n",
				i*f_bin, abs(ci), arg(ci)*M_180_PI, ci.real(), ci.imag(), harmonics[i]); 
		}   
		fclose(of);
	}

	if (F_wav || F_res)
	{	FILE* wavF = NULL;
		short* buf = NULL; // Buffer for raw sample data

		if (F_wav)
		{	if (strcmp(F_wav, "-") != 0)
			{	wavF = fopen(F_wav, "wb");
				if (wavF == NULL)
					die(41, "Failed to open %s for writing", F_wav);
			} else // stdout
			{	//_fsetmode(stdout, "b");
				wavF = stdout;
		}	}

		if (sweep)
		{	// Sweep mode
			fftw_real* sampbuf = new fftw_real[n_fft];

			if (F_wav)
			{	buf = new short[2*n_fft];
				if (n_rep)
					wavwriter(wavF, 2*n_fft * n_rep * fcount);
				else
					wavwriter(wavF, 0x1fffffdc);
			}

			// for each frequency
			for (size_t i = i_min; i <= i_max; ++i)
			{	if (fftbuf[i] == 0)
					continue; // skip this frequency
				// Design
				fprintf(stderr, "Freq: %f\t%f\n", f_bin * i, fftbuf[i]);
				for (size_t j = 0; j < n_fft; ++j)
					sampbuf[j] = cos(M_2PI * i * j / n_fft) * fftbuf[i];

				// write result
				if (F_res)
				{	FILE* of = fopen(F_res, resmode);
					if (of == NULL)
						die(41, "Failed to open %s for writing", F_res);
					const fftw_real* spe = sampbuf + n_fft;
					for (const fftw_real* sp = sampbuf; sp != spe; ++sp)
						fprintf(of, "%12g\n", *sp);
					fclose(of);
				}   

				if (execcmd)
					system(execcmd);

				// and quantize
				if (F_wav)
				{	quantize1(buf, sampbuf, n_fft);

					// Output
					size_t j = n_rep;
					do fwrite(buf, 2*n_fft * sizeof(short), 1, wavF);
					while (--j);

					// cleanup, well not in case of an exception...
				}
			} // for each frequency

			delete[] sampbuf;

		} else // if (sweep)
		{
			if (!stereo)
			{
				// ifft
				fftw_real* sampbuf = fftw_alloc_real(n_fft);
				fftw_plan plan = fftw_plan_r2r_1d(n_fft, fftbuf, sampbuf, FFTW_HC2R, FFTW_ESTIMATE);
				fftw_execute(plan);   // IFFT
				fftw_destroy_plan(plan);

				// normalize
				double fnorm = 0;
				fftw_real* sp = sampbuf;
				const fftw_real* const spe = sp + n_fft;
				for (; sp != spe; ++sp)
					if (fabs(*sp) > fnorm)
						fnorm = fabs(*sp);
				fnorm = 1/fnorm;
				for (sp = sampbuf; sp != spe; ++sp)
					*sp *= fnorm; 

				// write result
				if (F_res)
				{	FILE* of = fopen(F_res, resmode);
					if (of == NULL)
						die(41, "Failed to open %s for writing", F_res);
					for (sp = sampbuf; sp != spe; ++sp)
						fprintf(of, "%12g\n", *sp);
					fclose(of);
				}   

				// and quantize
				if (F_wav)
				{	buf = new short[2*n_fft];
					quantize1(buf, sampbuf, n_fft);
				}
				fftw_free(sampbuf); // no longer needed

			} else // stereo
			{
				// split channels
				fftw_real* const fftr = fftbuf + n_fft+1;
				for (size_t i = i_min; i <= i_max; ++i)
				{	if (harmonics[i] < 0)
						fftbuf[i] = fftbuf[n_fft-i] = 0;
				}

				// ifft
				fftw_real* sampbuf = fftw_alloc_real(2*n_fft);
				fftw_plan plan = fftw_plan_r2r_1d(n_fft, NULL, NULL, FFTW_HC2R, FFTW_ESTIMATE|FFTW_UNALIGNED);
				fftw_execute_r2r(plan, fftbuf, sampbuf);   // IFFT
				fftw_execute_r2r(plan, fftr, sampbuf+n_fft);
				fftw_destroy_plan(plan);

				// normalize
				{	double fnorm = 0;
					fftw_real* sp = sampbuf;
					const fftw_real* const spe = sp + 2*n_fft;
					for (; sp != spe; ++sp)
						if (fabs(*sp) > fnorm)
							fnorm = fabs(*sp);
					fnorm = 1/fnorm;
					for (sp = sampbuf; sp != spe; ++sp)
						*sp *= fnorm;
				}

				// write result
				if (F_res)
				{	FILE* of = fopen(F_res, resmode);
					if (of == NULL)
						die(41, "Failed to open %s for writing", F_res);
					const fftw_real* const spe = sampbuf + n_fft;
					for (fftw_real* sp = sampbuf; sp != spe; ++sp)
						fprintf(of, "%12g %12g %12g\n", sp[0]+sp[n_fft], sp[0], sp[n_fft]);
					fclose(of);
				}   

				// and quantize
				if (F_wav)
				{	buf = new short[2*n_fft];
					quantize2(buf, sampbuf, sampbuf + n_fft, n_fft);
				}
				fftw_free(sampbuf); // no longer needed
			} // if (stereo)

			if (F_wav)
			{	if (n_rep)
					wavwriter(wavF, 2*n_fft * n_rep);
				else
					wavwriter(wavF, 0x1fffffdc);
				// Output loop
				do fwrite(buf, 2*n_fft * sizeof(short), 1, wavF);
				while (--n_rep);
				// cleanup, well not in case of an exception...
			}

		} // if (sweep)

		delete[] buf;

		if (F_wav)
			fclose(wavF);
	} // if (F_wav || F_res)

 end:
	delete[] harmonics;
	fftw_free(fftbuf);
}

