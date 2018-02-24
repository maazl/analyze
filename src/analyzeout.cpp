#include "analyze.h"
#include "utils.h"
#include "mathx.h"
#include "pcmio.h"

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
#include <algorithm>
using namespace std;

#define M_2PI (2.*M_PI)
#define M_180_PI (180./M_PI)


AnalyzeOut::AnalyzeOut(const Config& cfg)
:	Cfg(cfg)
,	N2f((double)cfg.srate / cfg.N)
,	FminI((unsigned)ceil(cfg.fmin/N2f))
,	FmaxI((unsigned)fmin(floor(cfg.fmax/N2f), cfg.N/2))
,	OutLevel(fromdB(cfg.outgain))
,	LoopCount(!Cfg.zerooutfile ? Cfg.loops * Cfg.addloop : !Cfg.loops ? 0 : 2 * Cfg.loops * Cfg.addloop + Cfg.lpause + Cfg.discsamp/Cfg.N + 2)
,	PCMOut(cfg.floatsamp ? Format::F32 : cfg.swapbytes ? Format::I16_SWAP : Format::I16)
,	FCount(0)
{	FFTBuf.reset((cfg.stereo + 1) * (cfg.N + 1));
	OutBuf.reset(cfg.N * PCMOut.BytesPerSample);
}

bool AnalyzeOut::Setup()
{	//fprintf(stderr, "imin=%u imax=%u finc=%f flog=%f\n", FminI, FmaxI, cfg.f_inc, cfg.f_log);

	// generate coefficients
	unique_num_array<int> harmonics(Cfg.N + 2);
	FFTBuf.clear(); // all coefficients -> 0
	harmonics.clear();
	fftw_real maxamp = 0;
	int sign = 1;
	for (unsigned i = FminI; i <= FmaxI; ++i)
	{	if (i)
		{	// skip used harmonics
			for (unsigned j = 1; j <= Cfg.harmonic && i*j <= Cfg.N/2; ++j)
				if (harmonics[i*j])
					goto next_f; // continue in outer loop
			// lock harmonics
			for (unsigned j = 1; i*j <= Cfg.N/2; ++j)
				harmonics[i*j] = j * sign;
		}
		// calculate coefficients
		++FCount;
		FFTBuf[i] = pow(i, Cfg.scalepow);
		if (FFTBuf[i] > maxamp)
			maxamp = FFTBuf[i];
		if (!Cfg.sweep && i && i != Cfg.N/2)
		{	// apply random phase
			double phi = 2*M_PI * myrand();
			FFTBuf[Cfg.N-i] = FFTBuf[i] * sin(phi); // b[i]
			FFTBuf[i] *= cos(phi);                  // a[i]
		}
		// second channel
		if (sign < 0)
		{	FFTBuf[2*Cfg.N+1-i] = FFTBuf[Cfg.N-i];
			FFTBuf[Cfg.N+1+i] = FFTBuf[i];
		}
		// next frequency
		if (Cfg.stereo && !Cfg.sweep)
			sign = -sign;
		//fprintf(stderr, "f %i %i\n", i, (int)floor(i * f_log + f_inc));
		i = (unsigned)floor(i * Cfg.f_log + Cfg.f_inc - .5);
	 next_f:;
	}
	// normalize
	FFTBuf /= maxamp;

	// write design result
	if (Cfg.specfile)
	{	FILEguard of(Cfg.specfile, "w");
		fputs("#f\t|A|\targ A\tA real\tA imag\tharmon.\n", of);
		for (size_t i = 0; i <= Cfg.N/2; ++i)
		{	Complex ci(FFTBuf[i], i && i != Cfg.N/2 ? FFTBuf[Cfg.N-i] : 0);
			//          f     |A|  arg A ai  bi   harm
			fprintf(of, "%g\t%g\t%g\t%g\t%g\t%i\n",
				i*N2f, abs(ci), arg(ci)*M_180_PI, ci.real(), ci.imag(), harmonics[i]);
		}
	}

	if (!Cfg.outfile && !Cfg.reffile)
		return false; // nothing to do

	if (!Cfg.sweep)
	{	// Sweep mode needs late configuration for each frequency bin,
		// noise mode creates the PCM data here.
		if (!Cfg.stereo)
		{	// IFFT
			unique_fftw_arr<fftw_real> sampbuf(Cfg.N);
			fftwf_plan plan = fftwf_plan_r2r_1d(Cfg.N, FFTBuf.get(), sampbuf.get(), FFTW_HC2R, FFTW_ESTIMATE);
			fftwf_execute(plan);
			fftwf_destroy_plan(plan);

			Normalize(sampbuf);

			// write result
			if (Cfg.reffile)
			{	FILEguard of(Cfg.reffile, Cfg.refmode);
				for (auto v : sampbuf)
					fprintf(of, "%12g\n", v);
			}

			if (Cfg.outfile)
				PCMOut.convert(sampbuf, OutBuf);

		} else // stereo
		{	// split channels
			for (size_t i = FminI; i <= FmaxI; ++i)
				if (harmonics[i] < 0)
					FFTBuf[i] = FFTBuf[Cfg.N-i] = 0;

			// IFFT
			unique_fftw_arr<fftw_real> sampbuf(2*Cfg.N);
			fftwf_plan plan = fftwf_plan_r2r_1d(Cfg.N, NULL, NULL, FFTW_HC2R, FFTW_ESTIMATE|FFTW_UNALIGNED);
			fftwf_execute_r2r(plan, FFTBuf.get(), sampbuf.get());   // IFFT
			fftwf_execute_r2r(plan, FFTBuf.get() + Cfg.N + 1, sampbuf.get() + Cfg.N);
			fftwf_destroy_plan(plan);

			Normalize(sampbuf);

			// write result
			if (Cfg.reffile)
			{	FILEguard of(Cfg.reffile, Cfg.refmode);
				const fftw_real* const spe = sampbuf.begin() + Cfg.N;
				fputs("#l+r\tl\tr\n", of);
				for (fftw_real* sp = sampbuf.begin(); sp != spe; ++sp)
					fprintf(of, "%12g %12g %12g\n", sp[0]+sp[Cfg.N], sp[0], sp[Cfg.N]);
			}

			// and quantize
			if (Cfg.outfile)
				PCMOut.convert(sampbuf.slice(0, Cfg.N), sampbuf.slice(Cfg.N, Cfg.N), OutBuf);
		} // if (stereo)
	}

	if (Cfg.outfile)
		FOut = checkedopen(Cfg.outfile, "ab");

	return Cfg.outfile || (Cfg.reffile && Cfg.sweep);
}

void AnalyzeOut::Normalize(const unique_fftw_arr<fftw_real>& dst)
{
	double fnorm = 0;
	for (auto v : dst)
		if (fabs(v) > fnorm)
			fnorm = fabs(v);
	dst *= OutLevel / fnorm;
}

void AnalyzeOut::Run()
{
	if (Cfg.sweep)
	{	// Sweep mode
		unique_num_array<fftw_real> sampbuf(Cfg.N);

		if (FOut)
			PCMOut.WAVheader(FOut, Cfg.N * LoopCount * FCount, Cfg.srate);

		// for each frequency
		for (size_t i = FminI; i <= FmaxI && !termrq; ++i)
		{	double level = FFTBuf[i] * OutLevel;
			if (!level)
				continue; // skip this frequency
			// Design
			fprintf(stderr, "Freq: %f\t%f\n", N2f * i, level);
			double ff = M_2PI * i / Cfg.N;
			for (size_t j = 0; j < Cfg.N; ++j)
				sampbuf[j] = cos(ff * j) * level;

			// write result
			if (Cfg.reffile)
			{	FILEguard of(Cfg.reffile, Cfg.refmode);
				for (auto v : sampbuf)
					fprintf(of, "%12g\n", v);
			}

			/*if (execcmd)
				system(execcmd);*/

			// and play
			if (FOut)
			{	PCMOut.convert(sampbuf, OutBuf);
				PlayLoop();
			}
		} // for each frequency

	} else // if (sweep)
	{	// noise mode
		PCMOut.WAVheader(FOut, Cfg.N * LoopCount, Cfg.srate);
		PlayLoop();
	} // if (sweep)
}

void AnalyzeOut::PlayLoop()
{
	unsigned j = LoopCount;
	do
		fwrite2(OutBuf.begin(), OutBuf.size(), FOut);
	while (--j && !termrq);
}

