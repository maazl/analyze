#include "analyze.h"
#include "utils.h"
#include "mathx.h"
#include "pcmio.h"
#include "filereader.h"

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


fftw_real AnalyzeOut::MaxAbs(const unique_fftw_arr<fftw_real>& dst)
{	fftw_real ret = 0;
	for (auto v : dst)
	{	v = fabs(v);
		if (v > ret)
			ret = v;
	}
	return ret;
}

unsigned AnalyzeOut::CalcLoopCount(const Config& cfg)
{	if (!cfg.loops)
		return 0;
	unsigned count = cfg.loops * cfg.addloop;
	if (cfg.zerooutfile) // zero calibration: * 2 or 3
	{	if (cfg.normalize)
			count = (count * 2) + cfg.PauseLoops();
		else
			count = (count * 3) + 2 * cfg.PauseLoops();
	}
	count += (unsigned)ceil(cfg.predelay);
	if (!cfg.sweep && !cfg.sync)
		count += (cfg.discsamp + cfg.N) / cfg.N + 1;
	return count;
}

AnalyzeOut::AnalyzeOut(const Config& cfg)
:	Cfg(cfg)
,	N2f((double)cfg.srate / cfg.N)
,	MinQuant(Freq2EnergyQuantile(cfg.fmin))
,	MaxQuant(Freq2EnergyQuantile(cfg.fmax))
,	FminSmo(cfg.smooth ? EnergyQuantile2Freq(MinQuant + cfg.smooth * (MaxQuant-MinQuant) / 2.) : cfg.fmin)
,	FmaxSmo(cfg.smooth ? EnergyQuantile2Freq(MaxQuant - cfg.smooth * (MaxQuant-MinQuant) / 2.) : cfg.fmax)
,	FminI((unsigned)fmax(ceil((2*cfg.fmin - FminSmo)/N2f), 0.))
,	FmaxI(min((unsigned)floor((2*cfg.fmax - FmaxSmo)/N2f), cfg.N/2))
,	OutLevel(fromdB(cfg.outgain))
,	LoopCount(CalcLoopCount(cfg))
,	PCMOut(cfg.format, 1., (Cfg.outch == 0 && !(Cfg.stereo && (Cfg.sweep || Cfg.chirp))) - (Cfg.outch == 3))
{	Design.reset(cfg.N);
	Harmonics.reset(cfg.N/2 + 1);
	OutBuf.reset(cfg.N * PCMOut.BytesPerSample);
}

bool AnalyzeOut::Setup()
{	//fprintf(stderr, "imin=%u imax=%u finc=%f flog=%f\n", FminI, FmaxI, cfg.f_inc, cfg.f_log);
	FCount = 0;
	if (!Cfg.rspecfile)
		CreateDesign();
	else
		ReadDesign();

	unsigned secs = LoopCount;
	if (Cfg.sweep)
		secs *= FCount;
	if (Cfg.sweep || Cfg.chirp)
		secs <<= Cfg.stereo;
	secs = (unsigned)((uint64_t)secs * Cfg.N / Cfg.srate);
	fprintf(stderr, "Measurement time %u:%02u:%02u (loop count %u)\n",
		secs / 3600, secs / 60 % 60, secs % 60, LoopCount);

	// write design result
	if (Cfg.specfile)
	{	FILEguard of(Cfg.specfile, "w");
		fputs("#f\t|A|\targ A\tA real\tA imag\tharmon.\n", of);
		for (size_t i = 0; i <= Cfg.N/2; ++i)
		{	Complex ci(Design[i], i && i != Cfg.N/2 ? Design[Cfg.N-i] : 0);
			//          f     |A|  arg A ai  bi   harm
			fprintf(of, "%g\t%g\t%g\t%g\t%g\t%i\n",
				i*N2f, abs(ci), arg(ci)*M_180_PI, ci.real(), ci.imag(), Harmonics[i]);
		}
	}

	if (!Cfg.outfile && !Cfg.reffile)
		return false; // nothing to do

	// Sweep mode needs late configuration for each frequency bin,
	// noise mode creates the PCM data here. But for sync we need time domain data too.
	if (!Cfg.sweep || Cfg.sync)
	{	unique_fftw_arr<fftw_real> sampbuf;
		double norm;
		unique_fftwf_plan plan;

		CreateTimeDomain(sampbuf, norm, plan);

		CreatePCM(sampbuf);
	}

	if (Cfg.outfile)
		FOut = checkedopen(Cfg.outfile, "wb");

	return Cfg.outfile || (Cfg.reffile && Cfg.sweep);
}

double AnalyzeOut::Freq2EnergyQuantile(double f)
{	double xp = 1 + Cfg.scalepow / 2;
	if (xp == 0.)
		return log(f);
	return pow(f, xp) / xp;
}

double AnalyzeOut::EnergyQuantile2Freq(double q)
{	double xp = 1 + Cfg.scalepow / 2;
	if (xp == 0.)
		return exp(q);
	return pow(q * xp, 1 / xp);
}

void AnalyzeOut::CreateDesign()
{	Design.clear(); // all coefficients -> 0
	Harmonics.clear();

	// amplitude
	fftw_real maxamp = 0;
	int sign = 1;
	double sumamp = 0;
	for (unsigned i = FminI; i <= FmaxI; ++i)
	{	double f = i * N2f;
		if (i && !Cfg.sweep && !Cfg.chirp)
		{	// skip used harmonics
			for (unsigned j = 2; j < Cfg.harmonic && i*j <= Cfg.N/2; ++j)
				if (Harmonics[i*j])
					goto next_f; // continue in outer loop
			// lock harmonics
			for (unsigned j = 2; j < Cfg.harmonic && i*j <= Cfg.N/2; ++j)
				Harmonics[i*j] = j * sign;
		}
		Harmonics[i] = sign;
		++FCount;

		if (f < Cfg.fmin)
			f = 2 * Cfg.fmin - f;
		else if (f > Cfg.fmax)
			f = 2 * Cfg.fmax - f;
		Design[i] = pow(f, Cfg.scalepow / 2);
		//fprintf(stderr, "%i\t%g\t%g\t%g\t%g\t%g\n", i, i * N2f, f, Design[i], FminSmo, FmaxSmo);
		if (f < FminSmo)
			Design[i] *= .5 + .5 * sin(M_PI_2 * (i * N2f - Cfg.fmin) / (FminSmo - Cfg.fmin));
		else if (f > FmaxSmo)
			Design[i] *= .5 + .5 * sin(M_PI_2 * (i * N2f - Cfg.fmax) / (FmaxSmo - Cfg.fmax));
		sumamp += sqr(Design[i]);
		if (Design[i] > maxamp)
			maxamp = Design[i];

		// next frequency
		if (Cfg.stereo & !Cfg.sweep & !Cfg.chirp)
			sign = -sign;
		//fprintf(stderr, "f %i %i %i\n", i, (unsigned)floor(i * Cfg.f_log + Cfg.f_inc - .5), sign);
		i = (unsigned)floor(i * Cfg.f_log + Cfg.f_inc - .5);
	 next_f:;
	}

	// phase
	double phi = 0;
	double quantil = 0;
	for (unsigned i = max(FminI, 1U); i <= min(FmaxI, Cfg.N/2 - 1); ++i)
	{	quantil += sqr(Design[i]);
		if (!Cfg.chirp)
			phi = M_2PI * myrand();
		else
			phi += -M_2PI * quantil / sumamp * Cfg.chirp;
		Design[Cfg.N-i] = Design[i] * sin(phi); // b[i]
		Design[i] *= cos(phi);                  // a[i]
	}

	// normalize
	Design /= maxamp;
}

void AnalyzeOut::ReadDesign()
{	FileReader reader(Cfg.rspecfile);
	unique_array<double> data(6);
	for (unsigned bin = 0;; ++bin)
	{	if(!reader.ReadLine(data))
			die(30, "Design file does not have enough data for the current FFT size.\n"
				"Expected %u lines, EOF after %u lines.", Cfg.N/2 + 1, bin);
		if (fabs(data[0] - bin * N2f) > data[0] * 1E-5)
			die(30, "Design file data does not match the current FFT size or sampling frequency.\n"
				"Expected frequency: %f, found frequency %f.", bin * N2f, data[0]);
		fftw_real* design = Design.get();
		if (abs(Harmonics[bin]) == 1)
			++FCount;
		if ((Harmonics[bin] = (int)data[5]) < 0 && Cfg.stereo && !Cfg.sweep)
			design += Cfg.N + 1;
		if (bin)
			design[Cfg.N - bin] = data[4];
		if (bin == Cfg.N/2)
			return;
		design[bin] = data[3];
	}
}

void AnalyzeOut::CreateTimeDomain(unique_fftw_arr<fftw_real>& sampbuf, double& norm, unique_fftwf_plan& plan)
{
	if (!Cfg.stereo || Cfg.sweep || Cfg.chirp)
	{	sampbuf.reset(Cfg.N);
		plan = fftwf_plan_r2r_1d(Cfg.N, Design.get(), sampbuf.get(), FFTW_HC2R, FFTW_ESTIMATE);
		fftwf_execute(plan);
	} else // stereo
	{	sampbuf.reset(2*Cfg.N);
		unique_fftw_arr<fftw_real> design(Design.size());
		plan = fftwf_plan_r2r_1d(Cfg.N, design.get(), sampbuf.get(), FFTW_HC2R, FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// split channels
		design = Design;
		for (size_t i = FminI; i <= FmaxI; ++i)
			if (Harmonics[i] < 0)
			{	design[i] = 0;
				if (i && i < Cfg.N/2)
					design[Cfg.N-i] = 0;
			}
		fftwf_execute_r2r(plan, design.get(), sampbuf.get());
		design = Design;
		for (size_t i = FminI; i <= FmaxI; ++i)
			if (Harmonics[i] > 0)
			{	design[i] = 0;
				if (i && i < Cfg.N/2)
					design[Cfg.N-i] = 0;
			}
		fftwf_execute_r2r(plan, design.get(), sampbuf.get() + Cfg.N);
	}

	// normalize
	sampbuf *= OutLevel / (norm = MaxAbs(sampbuf));
}

void AnalyzeOut::CreatePCM(const unique_num_array<fftw_real>& sampbuf)
{
	if (sampbuf.size() == Cfg.N)
	{	if (Cfg.reffile)
		{	FILEguard of(Cfg.reffile, Cfg.refmode);
			for (auto v : sampbuf)
				fprintf(of, "%12g\n", v);
		}
		// and quantize
		if (Cfg.outfile)
			PCMOut.convert(sampbuf, OutBuf);
	} else
	{	if (Cfg.reffile)
		{	FILEguard of(Cfg.reffile, Cfg.refmode);
			const fftw_real* const spe = sampbuf.begin() + Cfg.N;
			fputs("#l+r\tl\tr\n", of);
			for (fftw_real* sp = sampbuf.begin(); sp != spe; ++sp)
				fprintf(of, "%12g %12g %12g\n", sp[0]+sp[Cfg.N], sp[0], sp[Cfg.N]);
		}
		// and quantize
		if (Cfg.outfile)
			PCMOut.convert(sampbuf.slice(0, Cfg.N), sampbuf.slice(Cfg.N, Cfg.N), OutBuf);
	}

	if (Cfg.outch == 2 && Cfg.outfile)
		NextChannel(); // Start with channel 2
}

void AnalyzeOut::Run()
{
	if (FOut)
	{	// WAV header
		if (!Cfg.nohdr)
		{	unsigned count = LoopCount;
			if (Cfg.sweep)
				count *= FCount;
			if (Cfg.sync)
				count += Cfg.sync + 1;
			PCMOut.WAVheader(FOut, Cfg.N * count, Cfg.srate);
		}

		// Sync preamble
		if (Cfg.sync)
		{	PlayNoise(Cfg.sync);
			// and one cycle silence as marker
			PlaySilence(1);
		}
	}

	if (Cfg.sweep)
		// Sweep mode
		DoSweep();
	else
	{	// noise or chirp mode
		PlayNoise(LoopCount);
		if (Cfg.stereo && Cfg.chirp)
		{	NextChannel();
			PlayNoise(LoopCount);
		}
	}

	if (FOut)
		fflush(FOut);
}

void AnalyzeOut::DoSweep()
{	unique_num_array<fftw_real> sampbuf(Cfg.N);

	// for each frequency
	for (size_t i = FminI; i <= FmaxI && !termrq; ++i)
	{	if (abs(Harmonics[i]) != 1)
			continue; // skip this frequency
		double level = sqrt(sqr(Design[i]) + sqr(Design[Cfg.N-i])) * OutLevel;
		// Design
		//fprintf(stderr, "Frequency: %f\t%f\n", N2f * i, level);
		sampbuf.clear();
		double ff = M_2PI * i / Cfg.N;
		for (size_t j = 0; j < Cfg.N; ++j)
			sampbuf[j] = sin(ff * j) * level;

		CreatePCM(sampbuf);

		if (FOut)
		{	PlayNoise(LoopCount);
			// 2nd channel
			if (Cfg.stereo)
			{ NextChannel();
				PlayNoise(LoopCount);
			}
		}
	} // for each frequency
}

void AnalyzeOut::PlayNoise(unsigned loopcount)
{	do
		fwrite2(OutBuf.begin(), OutBuf.size(), FOut);
	while (--loopcount && !termrq);
}

void AnalyzeOut::PlaySilence(unsigned loopcount)
{	loopcount *= OutBuf.size();
	char buf[1024] = {0};
	while (loopcount > 1024 && !termrq)
	{	fwrite2(buf, sizeof buf, FOut);
		loopcount -= 1024;
	}
	fwrite2(buf, loopcount, FOut);
}

void AnalyzeOut::NextChannel()
{	size_t shift = PCMOut.BytesPerSample >> 1;
	memmove(OutBuf.get() + shift, OutBuf.get(), OutBuf.size() - shift);
	memset(OutBuf.get(), 0, shift);
}
