#include "analyze.h"
#include "mathx.h"
#include "pca.h"
#include "utils.h"

#include <memory>
#include "filereader.h"

using namespace std;


void AnalyzeIn::AggEntry::add(double f, Complex I, Complex* U, unsigned harm, WeightFn wfn)
{	double w = wfn(abs(*U), abs(I), f);
	F += f * w;
	Iabs += abs(I) * w;
	Iarg += lIarg.Unwrap(arg(I)) * w;
	VE* val = Val;
	for (unsigned i = 0; i < harm; ++i, ++U, ++val)
	{	w = wfn(abs(*U), abs(I), f);
		val->Uabs += abs(*U) * w;
		val->Uarg += lUarg[i].Unwrap(arg(*U)) * w;
		Complex Z(*U / I);
		val->Zabs += abs(Z) * w;
		val->D += lZarg[i].Unwrap(arg(Z), f - lF) * w;
		val->Zarg += lZarg[i].Phase * w;
		val->W += w;
	}
	lF = f;
}

void AnalyzeIn::AggEntry::finish(unsigned harm)
{
	VE* val = Val;
	double w = val->W;
	F /= w;
	Iabs /= w;
	Iarg /= w;
	for (unsigned i = 0; i < harm; ++i, ++val)
	{	w = val->W;
		val->Uabs /= w;
		val->Uarg /= w;
		val->Zabs /= w;
		val->Zarg /= w;
		val->D    /= w * M_2PI;
	}
}


AnalyzeIn::FFTWorker::StoreRet AnalyzeIn::FFTWorker::StoreBin(const unsigned bin)
{
	// frequency
	double f = bin * Parent.N2f;
	if (f < Parent.Cfg.fmin)
		return BelowMin;
	if (f > Parent.Cfg.fmax)
		return AboveMax;

	Ch = Parent.SD.Harmonics[bin] < 0;
	FFTAgg& curagg = Agg[Ch];
	// retrieve coefficients
	Complex I(Parent.FFTBuffer[1][bin], bin && bin != Parent.Cfg.N / 2 ? Parent.FFTBuffer[1][Parent.Cfg.N - bin] : 0);
	Complex U[HA_MAX];
	// phase correction
	Complex phase(cos(Parent.LinPhase * f), sin(Parent.LinPhase * f));
	// calc U per harmonic
	curagg.Harm = 0;
	for (unsigned hb = bin; curagg.Harm < Parent.Cfg.harmonic && hb < Parent.Cfg.N / 2; hb += bin)
		U[curagg.Harm++] = Complex(Parent.FFTBuffer[0][hb], hb && hb != Parent.Cfg.N / 2 ? Parent.FFTBuffer[0][Parent.Cfg.N - hb] : 0) * phase;

	if (curagg.Bins == 0)
	{	// init
		curagg.next();
		curagg.NextF = f * (1 + Parent.Cfg.fbinsc) - Parent.N2f;
	}

	curagg.add(f, I, U, curagg.Harm, Parent.Cfg.weightfn);

	if (!curagg.Val[0].W)
		return Aggregated; // no weight so far => discard always

	++curagg.Bins;

	if (f < curagg.NextF)
		return Aggregated;

	curagg.finish(curagg.Harm);
	curagg.Bins = 0;
	//Zcache = polar(curagg->Zabs, curagg->Zarg);
	return Ready;
}

void AnalyzeIn::PrintHdr(FILE* dst)
{
  //     1   2    3      4    5      6    7      8       9       10      11     12
	fputs("#f\t|U|\targ U\t|I|\targ I\t|Z|\targ Z\tZ real\tZ imag\tweight\tdelay\tchan", dst);
	for (unsigned h = 1; h < Cfg.harmonic; ++h)
		fprintf(dst, "\t|Z%u|\targ Z%u\tZ%u real\tZ%u imag", h,h,h,h);
	fputc('\n', dst);
}

void AnalyzeIn::FFTWorker::PrintBin(FILE* dst) const
{
	const FFTAgg& curagg = Agg[Ch];
	Complex Z = polar(curagg.Val[0].Zabs, curagg.Val[0].Zarg);
	fprintf(dst, "%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%6i",
		// f      |Hl|                arg l                          |Hr|         arg r
		curagg.F, curagg.Val[0].Uabs, curagg.Val[0].Uarg * M_180_PI, curagg.Iabs, curagg.Iarg * M_180_PI,
		// |Z|              arg Z                          Z re      Z im
		curagg.Val[0].Zabs, curagg.Val[0].Zarg * M_180_PI, Z.real(), Z.imag(),
		// weight        delay            channel
		curagg.Val[0].W, curagg.Val[0].D, Ch);
	for (unsigned h = 1; h < Parent.Cfg.harmonic; ++h)
	{	const AggEntry::VE& val = curagg.Val[h];
		Z = polar(val.Zabs, val.Zarg);
		fprintf(dst, "\t%8g\t%8g\t%8g\t%8g",
		// |Hl|/|Hr| phil-phir           re        im
			val.Zabs, val.Zarg * M_180_PI, Z.real(), Z.imag());
	}
	fputc('\n', dst);
}


void AnalyzeIn::ReadColumn(const filecolumn& src, const unique_fftw_arr<fftw_real>& dst)
{	FILEguard in(src.file, "r");
	size_t count = dst.size();
	fftw_real* dp = dst.get();
	while (count--)
	{	ignore_result(fscanf(in, "#%*[^\n]")); // skip comments
		unsigned col = src.column;
		while (--col)
			ignore_result(fscanf(in, "%*s"));
		if (fscanf(in, "%g%*[^\n]", dp) != 1)
			die(27, "Failed to read column %u from data file.", src.column);
		++dp;
	}
}

void AnalyzeIn::CreateWindow(const unique_num_array<fftw_real>& dst, int type)
{	size_t len = dst.size() + 1;
	double sum = 0;
	fftw_real* win = dst.begin();
	for (size_t i = 1; i < len; i++)
	{	double w;
		switch (type)
		{default: // rectangle
			w = 1;
			break;
		 case 1: // Bartlett
			w = abs(i - len / 2.);
			break;
		 case 2: // Hanning
			w = .5 - .5 * cos(2 * M_PI * i / len);
			break;
		 case 3: // Hamming
			w = .54 - .46 * cos(2 * M_PI * i / len);
			break;
		 case 4: // Blackman
			w = .42 - .5 * cos(2 * M_PI * i / len) + .08 * cos(4 * M_PI * i / len);
			break;
		 case 5: // Blackman-Harris
			w = .35875 - .48829 * cos(2 * M_PI * i / len) + .14128 * cos(4 * M_PI * i / len) - .01168 * cos(6 * M_PI * i / len);
		}
		sum += *win++ = w;
	}
	dst *= len / sum;
}

void AnalyzeIn::ApplyCalibration()
{
	fftw_real* a1 = FFTBuffer[0].get();
	fftw_real* a2 = FFTBuffer[1].get();
	fftw_real* b1 = a1 + Cfg.N;
	fftw_real* b2 = a2 + Cfg.N;
	for (int bin = 0; a1 < b1; ++bin, ++a1, ++a2, --b1, --b2)
	{	// calc
		double f = bin * N2f;
		Complex U(bin != Cfg.N / 2 ? *a1 : 0, bin ? *b1 : 0);
		Complex I(bin != Cfg.N / 2 ? *a2 : 0, bin ? *b2 : 0);
		//fprintf(stderr, "%g, *a1=%g, *a2=%g, *b1=%g, *b2=%g\n", f, *a1, *a2, *b1, *b2);
		// calibration
		if (Cfg.gaininfile && !Cfg.gainoutfile)
			U /= Gain[bin];
		if (Cfg.gainoutfile)
		{	double weight = sqrt(Cfg.weightfn(abs(U), abs(I), f));
			WSums[bin] += weight;
			U /= Gain[bin];
			GainD[bin] += I / U * weight;
		}
		if (Cfg.zeroinfile)
		{	const MatrixC<2,2>& z = Zero[bin];
			Complex t = U;
			U = (U * z[0][0] + I * z[1][0]);
			I = (t * z[0][1] + I * z[1][1]);
		}
		if (Cfg.zerooutfile)
		{	VectorC<2>& z = ZeroD[bin][!ZeroPart2];
			z[0] += U;
			z[1] += I;
		}
		// store data
		if (bin != Cfg.N / 2)
		{	*a1 = U.real();
			*a2 = I.real();
		}
		if (bin)
		{	*b1 = U.imag();
			*b2 = I.imag();
		}
	}
}

AnalyzeIn::AnalyzeIn(const Config& cfg, const SetupData& sd)
	// setup input
:	Cfg(cfg)
,	SD(sd)
,	N2f((double)Cfg.srate / Cfg.N)
,	PCMIn(Cfg.floatsamp ? Format::F32 : Cfg.swapbytes ? Format::I16_SWAP : Format::I16, Cfg.diffmode, &Cfg.gainadj)
,	LinPhase(Cfg.linphase * M_2PI)
{	// allocate buffers
	InBufferTmp.reset(PCMIn.BytesPerSample * Cfg.N);
	Window.reset(Cfg.N);
	InBuffer[0].reset(Cfg.N);
	InBuffer[1].reset(Cfg.N);
	if (Cfg.mxy)
	{	// reserve space for integrals and differentials too
		IntBuffer[0].reset(Cfg.N);
		IntBuffer[1].reset(Cfg.N);
		DiffBuffer[0].reset(Cfg.N);
		DiffBuffer[1].reset(Cfg.N);
		FFTBuffer[0].reset(3 * Cfg.N + 1);
		FFTBuffer[1].reset(3 * Cfg.N + 1);
	} else
	{	FFTBuffer[0].reset(Cfg.N + 1);
		FFTBuffer[1].reset(Cfg.N + 1);
	}
	if (Cfg.crosscorr | Cfg.sync)
	{	CCBuffer1.reset(Cfg.N);
		CCBuffer2.reset(Cfg.N);
	}
	WSums.reset(Cfg.N / 2 + 1);
	WSums.clear();
	Gain.reset(Cfg.N / 2 + 1);
	GainD.reset(Cfg.N / 2 + 1);
	Zero.reset(Cfg.N / 2 + 1);
	ZeroD.reset(Cfg.N / 2 + 1);
	// create plan
	// fftw_real in[N], tout[N], power_spectrum[N/2+1];
	P = fftwf_plan_r2r_1d(Cfg.N, InBuffer[0].get(), FFTBuffer[0].get(), FFTW_R2HC, FFTW_ESTIMATE);
	PI = fftwf_plan_r2r_1d(Cfg.N, InBuffer[0].get(), FFTBuffer[0].get(), FFTW_HC2R, FFTW_ESTIMATE);
}

bool AnalyzeIn::Setup()
{	// adjust fmax
	double fmax = min(Cfg.fmax, (double)Cfg.srate/2);

	CreateWindow(Window, Cfg.winfn);
	// write window data
	if (Cfg.windowfile)
	{	FILEguard out(Cfg.windowfile, "wt");
		for (fftw_real val : Window)
			fprintf(out, "%g\n", val);
	}

	if (Cfg.overwrt[0].file)
	{	OvrBuffer[0].reset(Cfg.N);
		ReadColumn(Cfg.overwrt[0], OvrBuffer[0]);
	}
	if (Cfg.overwrt[1].file)
	{	OvrBuffer[1].reset(Cfg.N);
		ReadColumn(Cfg.overwrt[1], OvrBuffer[1]);
	}

	GainD.clear();
	ZeroD.clear();
	// read gain file
	if (Cfg.gaininfile)
	{	//PolarFileInterpolation ip(Cfg.gaininfile, 1);
		FileInterpolation ip(Cfg.gaininfile, 2);
		unsigned i = 0;
		for (auto& g : Gain)
		{	auto& row = ip.Get(i++ * N2f);
			g = Complex(row[1], row[2]);
		}
	}
	// read zero file
	if (Cfg.zeroinfile)
	{	//PolarFileInterpolation ip(Cfg.zeroinfile, 4);
		FileInterpolation ip(Cfg.zeroinfile, 8);
		unsigned i = 0;
		for (auto& z : Zero)
		{	auto& row = ip.Get(i++ * N2f);
			z[0][0] = Complex(row[1], row[2]);
			z[0][1] = Complex(row[3], row[4]);
			z[1][0] = Complex(row[5], row[6]);
			z[1][1] = Complex(row[7], row[8]);

			if (Cfg.normalize)
			{	Complex s = 1. / (z[0][0] + z[1][0]);
				z[0][0] *= s;
				z[1][0] *= s;
				s = 1. / (z[0][1] + z[1][1]);
				z[0][1] *= s;
				z[1][1] *= s;
			}
			/*fprintf(stderr, "Zero %f\t%p\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (i-1) * N2f, &row,
				z[0][0].real(), z[0][0].imag(), z[0][1].real(), z[0][1].imag(), z[1][0].real(), z[1][0].imag(), z[1][1].real(), z[1][1].imag());*/
			z = inverse(z); // invert matrix to save time during analysis
			//fprintf(stderr, "ZeroI\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", z[0][0], z[0][1], z[1][0], z[1][1]);
		}
	}

	// Calculate autocorrelation for SNR calculations of cross correlations
	if (Cfg.crosscorr | Cfg.sync)
	{	SNRAC = abs(ExecuteCrossCorrelation(SD.Design, SD.Design));
		fprintf(stderr, "Autocorrelation level of reference: %g\n", SNRAC);
	}

	if (Cfg.infile)
		FIn = checkedopen(Cfg.infile, "rb");

	return Cfg.infile || (Cfg.overwrt[0].file && Cfg.overwrt[1].file);
}

void AnalyzeIn::Run()
{
	if (Cfg.infile)
		// discard first samples
		PCMIn.discard(FIn, Cfg.discsamp, InBufferTmp);

	NeedSync = Cfg.sync || Cfg.predelay;
	NextSamples = Cfg.N;
	SyncPhaseVector = 0;
	SyncWeight = 0;

	unique_ptr<SweepWorker> sweep;
	if (Cfg.sweep)
		sweep.reset(new SweepWorker(*this));

 restart_zero:
	// operation loop
	unsigned loop = Cfg.loops * Cfg.addloop;
	if (Cfg.sweep)
		loop *= SD.FCount;
	unsigned addloops = 0;
	unsigned block = 0;
	do
	{
		if (FIn)
		{resync:
			fprintf(stderr, "loop block = %u, nsync = %i, nextsamp = %u\n", block, NeedSync, NextSamples);
			if (sweep && !addloops && !NeedSync && !block)
				// frequency setup delay for sweep mode
				PCMIn.discard(FIn, Cfg.N * (unsigned)ceil(Cfg.predelay), InBufferTmp);

			// Read NextSamples and keep Cfg.N - NextSamples if > 0.
			const unsigned overlap = (Cfg.N - NextSamples);
			// read data
			const auto& in = InBufferTmp.slice(0, InBufferTmp.size() - overlap * PCMIn.BytesPerSample);
			fread2(in.get(), in.size(), FIn);
			// write raw data
			if (Cfg.rawfile)
				PCMIn.ASCIIdump(FILEguard(Cfg.rawfile, "wt"), in);

			// Keep some samples?
			if (overlap)
			{	InBuffer[0].slice(0, overlap) = InBuffer[0].slice(NextSamples, overlap);
				InBuffer[1].slice(0, overlap) = InBuffer[1].slice(NextSamples, overlap);
			}
			// reset min/max
			PCMIn.reset();
			// convert data
			PCMIn.convert(InBuffer[Cfg.swapch].slice(overlap, NextSamples), InBuffer[!Cfg.swapch].slice(overlap, NextSamples), in, !NeedSync && (addloops || (Cfg.incremental && block)));
			// write raw status
			if (PCMIn.Limits[0].Min == PCMIn.Limits[0].Max && PCMIn.Limits[1].Min == PCMIn.Limits[1].Max)
				fprintf(stderr, "Input data is zero, check mute switch!\n");
			if (Cfg.minmax)
				fprintf(stderr, "\nmin:\t%g\t%g\nmax:\t%g\t%g\n",
						PCMIn.Limits[0].Min, PCMIn.Limits[1].Min, PCMIn.Limits[0].Max, PCMIn.Limits[1].Max);

			// synchronize?
			if (NeedSync)
			{	DoSync();
				if (termrq)
					break;
				goto resync;
			}
			NextSamples = Cfg.N;
		}

		if (addloops || (Cfg.incremental && block))
		{	if (OvrBuffer[0])
				InBuffer[0] += OvrBuffer[0];
			if (OvrBuffer[1])
				InBuffer[1] += OvrBuffer[1];
		} else
		{	if (OvrBuffer[0])
				InBuffer[0] = OvrBuffer[0];
			if (OvrBuffer[1])
				InBuffer[1] = OvrBuffer[1];
		}

		//fprintf(stderr, "AL %i/%i\n", addloops, Cfg.addloop);
		if (++addloops < Cfg.addloop)
			continue; // add more data
		addloops = 0;

		// write source data
		if (Cfg.srcfile)
		{	FILEguard out(Cfg.srcfile, "wt");
			size_t len = InBuffer[0].size();
			const fftw_real* sp1 = InBuffer[0].get();
			const fftw_real* sp2 = InBuffer[1].get();
			while (len--)
				fprintf(out, "%g\t%g\n", *sp1++, *sp2++);
		}

		if (Cfg.winfn)
		{	InBuffer[0] *= Window;
			InBuffer[1] *= Window;
		}

		// Do the main processing
		if (sweep)
			sweep->DoSweep(block);
		else if (Cfg.mfft & Cfg.mpca)
			DoFFTPCA();
		else if (Cfg.mpca)
			DoPCA();
		else if (Cfg.mfft)
			DoFFT();
		else if (Cfg.mxy)
			DoXY();

		Cfg.plot.execute();

		if (++block == Cfg.loops)
			block = 0;;
	} while (--loop && !termrq);

	if (Cfg.gainoutfile)
	{	FILEguard fz(Cfg.gainoutfile, "wb");
		fputs("#f\treal\timag\tabs\targ\n", fz);
		unsigned i = 0;
		double* wp = WSums.begin();
		for (Complex g : GainD)
		{	g /= *wp++;  // scale 2 average
			fprintf(fz, "%f\t%g\t%g\t%g\t%g\n", i++ * N2f, g.real(), g.imag(), abs(g), arg(g) * M_180_PI);
		}
	}
	if (Cfg.zerooutfile)
	{	if (!ZeroPart2) // part 1
		{	ZeroPart2 = true;
			fputs("Zeromode calibration part one is now complete.\n"
				"Part 2 will start at the end of the contdown.\n\7", stderr);
			for (unsigned loop = Cfg.PauseLoops(); loop; --loop)
			{	fprintf(stderr, "\r%u ", loop);
				PCMIn.discard(FIn, Cfg.N, InBufferTmp);
			}
			puts("\nNow at part 2.");
			goto restart_zero;
		} else // part 2
		{	FILEguard fz(Cfg.zerooutfile, "wb");
			fputs("#f\tU->U re\tU->U im\tU->I re\tU->I im\tI->U re\tI->U im\tI->I re\tI->I im\t" // 0..8
				"|U->U|\targ U->U\t|U->I|\targ U->I\t|I->U|\targ I->U\t|I->I|\targ I->I\t" // 9..16
				"|Uo|\targ Uo\t|U0|\targ U0\t|Io|\targ Io\t|I0|\targ I0\n", fz); // 17..24
			unsigned i = 0;
			for (const auto& z : ZeroD)
			{	auto zn = z;
				if (Cfg.normalize)
				{	Complex s = 1. / (zn[0][0] + zn[1][0]);
					zn[0][0] *= s;
					zn[1][0] *= s;
					s = 1. / (zn[0][1] + zn[1][1]);
					zn[0][1] *= s;
					zn[1][1] *= s;
				}
				// scale to fit det z == 1
				Complex sca = 1. / sqrt(det(zn));
				if (((zn[0][0] + zn[1][1]) * sca).real() < 0)
					sca = -sca;
				auto z_ = zn * sca;
				fprintf(fz, "%f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					i++ * N2f,
					z_[0][0].real(), z_[0][0].imag(), z_[0][1].real(), z_[0][1].imag(),
					z_[1][0].real(), z_[1][0].imag(), z_[1][1].real(), z_[1][1].imag(),
					abs(z_[0][0]), arg(z_[0][0]) * M_180_PI, abs(z_[0][1]), arg(z_[0][1]) * M_180_PI,
					abs(z_[1][0]), arg(z_[1][0]) * M_180_PI, abs(z_[1][1]), arg(z_[1][1]) * M_180_PI,
					abs(z[0][0]), arg(z[0][0]) * M_180_PI, abs(z[0][1]), arg(z[0][1]) * M_180_PI,
					abs(z[1][0]), arg(z[1][0]) * M_180_PI, abs(z[1][1]), arg(z[1][1]) * M_180_PI);
			}
		}
	}

	fputs("completed.", stderr);

	Cfg.post.execute();

	// read until EOF
	if (Cfg.disctrail && !termrq)
		while (fread(InBufferTmp.get(), sizeof InBufferTmp[0], InBufferTmp.size(), FIn) > 0);
}

void AnalyzeIn::DoSync()
{
	if (!Cfg.sync)
	{haveSync:
		NeedSync = false;
		NextSamples = Cfg.N;
		return;
	}

	// select channels
	fftw_real* in;
	switch (Cfg.syncch)
	{default:
		InBuffer[0].slice(Cfg.N - NextSamples, NextSamples) += InBuffer[1].slice(Cfg.N - NextSamples, NextSamples);
	 case 1:
		in = InBuffer[0].get();
		break;
	 case 2:
		in = InBuffer[1].get();
	}
	// forward transformation
	fftwf_execute_r2r(P, in, FFTBuffer[0].get());

	// append constant zero to FFT result to simplify analysis
	FFTBuffer[0][Cfg.N] = 0;

	// cross correlate with design function
	Complex ccv = ExecuteCrossCorrelation(FFTBuffer[0], SD.Design) / SNRAC;
	double level = abs(ccv);
	if (SyncWeight > 0 && level < abs(SyncPhaseVector) / SyncWeight * Cfg.syncend)
	{	fprintf(stderr, "End of sync signal after %f cycles.\n", SyncWeight);
		goto haveSync;
	} else if (!isfinite(level) || level < Cfg.synclevel)
	{	if (SyncWeight < Cfg.sync * -1.5)
			die(29, "Failed to synchronize after %f cycles.\n", 1-SyncWeight);
		else
		{	fprintf(stderr, "No sync signal: cross correlation level %f too low.\n", level);
			--SyncWeight;
			NextSamples = Cfg.N; // try next N samples
			return;
		}
	}
	// Now abs(ccv) is the relative sync level of the last N samples
	// and arg(ccv) is the phase shift relative to the design function.

	// Aggregate cross correlation result
	// Calculate delay, i.e. NextSamples
	double shift = arg(SyncPhaseVector + ccv) / M_2PI;
	NextSamples = (unsigned)(Cfg.N * (shift + (shift <= 0)) + .5);
	if (NextSamples == 0)
		NextSamples = Cfg.N;
	// calculate weight
	double rounded_shift = (double)NextSamples / Cfg.N;
	SyncPhaseVector += ccv * rounded_shift; // No. of new samples added this cycle is the weight.
	if (SyncWeight < 0)
		SyncWeight = 0;
	SyncWeight += rounded_shift;

	fprintf(stderr, "Synchronized with time shift %f samples = %g ms, correlation level %f\n",
		shift * Cfg.N, 1E3 * shift / N2f, abs(SyncPhaseVector) / SyncWeight);

	// adjust SPV according to NextSamples
	// This leaves only the residuals of less than one sample in the phase of SPV.
	SyncPhaseVector *= polar(1., -M_2PI * rounded_shift);
	//fprintf(stderr, "SPV %f + %g i\n", SyncPhaseVector.real(), SyncPhaseVector.imag());
	fprintf(stderr, "request %u samples\n", NextSamples);
}

void AnalyzeIn::DoPCA()
{
	PCA<5> pca;
	double data[6];
	fftw_real* U = InBuffer[0].get() + 1;
	fftw_real* I = InBuffer[1].get() + 1;
	data[2] = 1;   // offset
	data[3] = 0;   // integral
	data[4] = 0;   // linear
	data[5] = 0;   // differential
	unsigned i = Cfg.N - 3;
	do
	{	data[0] = U[0] + U[1];
		data[1] = I[0] + I[1];
		data[3] += I[-1] + I[0];
		data[5] = I[-1] + I[0] - I[1] - I[2];
		pca.Store((double (&)[5])data, Cfg.weightfn(data[0], data[1], i));
		data[4]++;
		U += 2;
		I += 2;
		i -= 2;
	} while (i > 0);

	VectorD<4> res = pca.Result() * Cfg.rref;
	fprintf(stderr, "\nPCA: %12g %12g %12g %12g %12g %12g\n",
		res[0], res[1], 2. / Cfg.srate / res[2], res[3], 1. / 2 * Cfg.srate * res[2] / M_2PI / res[0], Cfg.srate * res[4] / 2);
}

void AnalyzeIn::DoFFT()
{
	ExecuteFFT();

	FFTWorker fft(*this);

	// calc some sums
	double wsum = 0;
	int nsum = 0;
	double Rsum = 0;
	double R2sum = 0;
	double L2sum = 0;
	//double LCsum = 0; ==> -2 * wsum
	double C2sum = 0;
	double Lsum = 0;
	double Csum = 0;
	double d2sum = 0;
	//double RWsum = 0;
	// write data
	FILEguard tout = NULL;
	if (Cfg.datafile)
	{	tout = checkedopen(Cfg.datafile, "wt");
		PrintHdr(tout);
	}

	for (size_t len = 0; len <= Cfg.N / 2; ++len)
	{	// do calculations and aggregations
		switch (fft.StoreBin(len))
		{default:
			//case FFTbin::Aggregated:
			continue;
		 case FFTWorker::Ready:
			// write
			if (tout)
				fft.PrintBin(tout);
		}

		const AggEntry& calc = fft.ret();
		if (calc.F < Cfg.famin || calc.F > Cfg.famax)
			continue;
		const Complex z = polar(calc.Val[0].Zabs, calc.Val[0].Zarg);
		//fprintf(stderr, "FW %12g %12g\n", calc.f(), calc.W());
		// average
		++nsum;
		wsum += calc.Val[0].W;
		// resistivity
		Rsum += calc.Val[0].W * z.real();
		R2sum += calc.Val[0].W * sqr(z.real());
		// L & C
		L2sum += calc.Val[0].W * sqr(calc.F);
		// LCsum += weight; == wsum
		C2sum += calc.Val[0].W / sqr(calc.F);
		Lsum += calc.Val[0].W * z.imag() * calc.F;
		Csum += calc.Val[0].W * z.imag() / calc.F;
		d2sum = calc.Val[0].W * sqr(z.imag());
	}

	// calculate summary
	double R = Rsum / wsum;
	double RE = sqrt((R2sum / wsum - sqr(R)) / (nsum - 1));
	double sLC = 1 / (sqr(wsum) - C2sum * L2sum);
	double L = (C2sum * Lsum - Csum * wsum) * sLC;
	double C = (Csum * L2sum - wsum * Lsum) * sLC;
	double LE = sqrt((C2sum * d2sum - sqr(Csum) - (d2sum * sqr(wsum) - 2 * Csum * wsum * Lsum + C2sum * sqr(Lsum)) / L2sum) * sLC / (nsum - 2));
	double CE = sqrt((L2sum * d2sum - sqr(Lsum) - (d2sum * sqr(wsum) - 2 * Csum * wsum * Lsum + L2sum * sqr(Csum)) / C2sum) * sLC / (nsum - 2));

	// write summary
	// fprintf(stderr, "\n%6i %12g %12g %12g %12g %12g %12g %12g %12g\n", nsum, wsum, Rsum, R2sum, Csum, C2sum, Lsum, L2sum, sLC);
	fprintf(stderr,
		"\nreal: ESR [Ω]\t%12g ± %g\n"
		"imag.: ESC [µF]\t%12g ± %g\n"
		"imag.: ESL [µH]\t%12g ± %g\n",
		Cfg.rref * R, Cfg.rref * RE, 1E6 / (Cfg.rref * C * M_2PI), 1E6 / (CE * Cfg.rref * M_2PI),
		1E6 * Cfg.rref * L / M_2PI, 1E6 * Cfg.rref * LE / M_2PI);
	if (Cfg.crosscorr)
		fprintf(stderr, "delay [ms]\t%12g\n", 1E3 / M_2PI * LinPhase);
}

void AnalyzeIn::DoFFTPCA()
{	// FFT
	ExecuteFFT();

	FFTWorker fft(*this);

	// write data
	FILEguard tout = NULL;
	if (Cfg.datafile)
	{	tout = checkedopen(Cfg.datafile, "wt");
		PrintHdr(tout);
	}

	PCA<2> pcaRe;
	PCA<3> pcaIm;
	double PCAdataRe[2];
	double PCAdataIm[3];
	// some values are const
	PCAdataRe[1] = 1;
	//PCAdataIm[1] = 1;

	// 1st line
	for (size_t len = 0; len <= Cfg.N / 2; ++len)
	{	// do calculations and aggregations
		switch (fft.StoreBin(len))
		{default:
			//case FFTbin::Aggregated:
			//case FFTbin::Skip:
			continue;
		 case FFTWorker::Ready:
			// write
			if (tout)
				fft.PrintBin(tout);
		}

		const AggEntry& calc = fft.ret();
		if (calc.F < Cfg.famin || calc.F > Cfg.famax)
			continue;
		const Complex z = polar(calc.Val[0].Zabs, calc.Val[0].Zarg);
		// component analysis
		PCAdataRe[0] = z.real();
		//PCAdataRe[2] = 1/af;
		//PCAdataRe[3] = f;
		PCAdataIm[0] = z.imag(); // fit imaginary part in conductivity
		PCAdataIm[1] = 1. / calc.F;
		PCAdataIm[2] = calc.F;
		//PCAdataIm[3] = 1/af;
		//fprintf(stderr, "Re: %12g %12g %12g %12g\n", calc.f(), PCAdataRe[0], PCAdataRe[1], calc.W());

		// add values
		pcaRe.Store(PCAdataRe, calc.Val[0].W);
		pcaIm.Store(PCAdataIm, calc.Val[0].W);
	}

	// calculate summary
	VectorD<1> resRe = pcaRe.Result();
	VectorD<2> resIm = pcaIm.Result();

	//fprintf(stderr, "resRe: %12g %12g %12g\n", resRe[0], resRe[1], resRe[2]);
	//fprintf(stderr, "resIm: %12g %12g %12g %12g\n", resIm[0], resIm[1], resIm[2], resIm[3]);

	double R0 = Cfg.rref * resRe[0];
	//double R1_f = rref*resRe[1];
	double C0 = -1 / M_2PI / resIm[0] / Cfg.rref;
	double L0 = resIm[1] / M_2PI * Cfg.rref;
	//double C1f = -resIm[0] / M_2PI / rref;

	// write summary
	fprintf(stderr, "\nreal: R [Ohm]\t%12g\n"
	//       "      R/f [Ohm/Hz]\t%12g\t%12g @100Hz\n"
	    "imag.: C [µF]\t%12g\n"
	//  "      Cf [F Hz]   \t%12g\t%12g @100Hz\n"
	    "imag.: L [µH]\t%12g\n"
	//             , R0, R1_f, R1_f/100, C0, C1f, C1f*100);
	, R0, C0 * 1E6, L0 * 1E6);
	if (Cfg.crosscorr)
		fprintf(stderr, "delay [ms]\t%12g\n", 1E3 / M_2PI * LinPhase);
}

AnalyzeIn::SweepWorker::SweepWorker(AnalyzeIn& parent)
:	Parent(parent)
,	SweepFrequency(0)
,	Channel(0)
,	LastFrequency(0)
{	if (Parent.Cfg.datafile)
	{	Out = checkedopen(Parent.Cfg.datafile, "w");
		Parent.PrintHdr(Out);
	}
}

void AnalyzeIn::SweepWorker::DoSweep(unsigned block)
{	// find next frequency
	if (!SweepFrequency)
		SweepFrequency = (unsigned)ceil(Parent.Cfg.fmin/Parent.N2f); // start at fmin
	else if (block == 0 && ++Channel > Parent.Cfg.stereo)
	{	Channel = 0;
		LastFrequency = SweepFrequency;
		while (abs(Parent.SD.Harmonics[++SweepFrequency]) != 1)
			assert(SweepFrequency * Parent.N2f <= Parent.Cfg.fmax); // must not overflow and must not exceed fmax
	}
	fprintf(stderr, "Frequency: %f, bin: %u, channel: %u\n", Parent.N2f * SweepFrequency, SweepFrequency, Channel);

	// Analyze frequency including harmonics
	const double fs = -M_2PI * SweepFrequency / Parent.Cfg.N;
	const fftw_real* up = Parent.InBuffer[0].get();
	const fftw_real* ip = Parent.InBuffer[1].get();
	for (unsigned i = 0; i < Parent.Cfg.N; ++i)
	{	// calculate sin/cos sums
		auto phasevec = polar(1., fs * i);
		Ana[0].store(phasevec, *ip++);
		double val = *up++;
		Ana[1].store(phasevec, val);
		for (unsigned j = 2; j <= Parent.Cfg.harmonic; ++j)
			Ana[j].store(polar(1., fs * i * j), val);
	}
	// post process
	for (unsigned j = 0; j <= Parent.Cfg.harmonic; ++j)
		Ana[j].finish(Parent.Cfg.N);

	if (Out)
		Print();
}

void AnalyzeIn::SweepWorker::Print()
{
	double f = SweepFrequency * Parent.N2f;
	double U = abs(Ana[1].Xi);
	double I = abs(Ana[0].Xi);
	Complex Z = Ana[1].Xi / Ana[0].Xi;
	GroupDelay& delay = Delay[Channel];
	fprintf(Out, "%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%6i",
		// f |Hl| phil                   |Hr| phir
		f, U, arg(Ana[1].Xi) * M_180_PI, I, arg(Ana[0].Xi) * M_180_PI,
		// |Hl|/|Hr| phil-phir          re        im
		abs(Z), delay.Phase * M_180_PI, Z.real(), Z.imag(),
		// weight                     delay           channel
		Parent.Cfg.weightfn(U, I, f), delay.Unwrap(arg(Z), (SweepFrequency - LastFrequency) * Parent.N2f), Channel);
	Unwrapper* phase = Phase[Channel];
	for (unsigned h = 2; h <= Parent.Cfg.harmonic; ++h, ++phase)
	{	Z = Ana[h].Xi / Ana[0].Xi;
		fprintf(Out, "\t%8g\t%8g\t%8g\t%8g",
			// |Hl|/|Hr| phil-phir                    re        im
			abs(Z), phase->Unwrap(arg(Z)) * M_180_PI, Z.real(), Z.imag());
	}
	fputc('\n', Out);
	fflush(Out);
}

void AnalyzeIn::DoXY()
{	// we need to do an FFT here, at least for the calibration
	ExecuteFFT();

	// U(f)
	fftw_real* a1 = FFTBuffer[0].get();
	fftw_real* a2 = FFTBuffer[1].get();
	fftw_real* b1 = a1 + Cfg.N;
	fftw_real* b2 = a2 + Cfg.N;
	*b1 = 0; // well, somewhat easier this way
	*b2 = 0;
	for (int len = 0; a1 < b1; len += Cfg.harmonic, a1 += Cfg.harmonic, a2 += Cfg.harmonic, b1 -= Cfg.harmonic, b2 -= Cfg.harmonic)
	{  // calc
		double f = len * N2f;
		Complex U(*a1, *b1);
		Complex I(*a2, *b2);
		// The integrals and differentials are calculated in the frequency domain.
		// This is at the expence of approximate a factor 2 computing time, since we have to do
		// one forward and one backward transformation for the zero compensation anyway.
		// The advantage is that we do not have to deal with the 1/2 time slot linear phase shift
		// of the numeric integration/differentiation in the time domain.
		Complex di(0, f); // differential operator
		// store integrals
		Complex C = U / di;
		a1[Cfg.N] = C.real();
		b1[Cfg.N] = C.imag();
		C = I / di;
		a2[Cfg.N] = C.real();
		b2[Cfg.N] = C.imag();
		// store differentials
		C = U * di;
		a1[2 * Cfg.N] = C.real();
		b1[2 * Cfg.N] = C.imag();
		C = I * di;
		a2[2 * Cfg.N] = C.real();
		b2[2 * Cfg.N] = C.imag();
	}
	// clear DC and Nyquist frequencies of integral and differential, since they do not allow 90° phase shift.
	FFTBuffer[0][Cfg.N] = FFTBuffer[1][Cfg.N] = 0;
	FFTBuffer[0][2 * Cfg.N] = FFTBuffer[1][2 * Cfg.N] = 0;
	FFTBuffer[0][Cfg.N + Cfg.N / 2] = FFTBuffer[1][Cfg.N + Cfg.N / 2] = 0;
	FFTBuffer[0][2 * Cfg.N + Cfg.N / 2] = FFTBuffer[1][2 * Cfg.N + Cfg.N / 2] = 0;

	// now do the inverse transform to get the corrected data back.
	fftwf_execute_r2r(PI, FFTBuffer[0].get(), InBuffer[0].get());
	fftwf_execute_r2r(PI, FFTBuffer[0].get() + Cfg.N, IntBuffer[0].get());
	fftwf_execute_r2r(PI, FFTBuffer[0].get() + 2 * Cfg.N, DiffBuffer[0].get());
	fftwf_execute_r2r(PI, FFTBuffer[1].get(), InBuffer[1].get());
	fftwf_execute_r2r(PI, FFTBuffer[1].get() + Cfg.N, IntBuffer[1].get());
	fftwf_execute_r2r(PI, FFTBuffer[1].get() + 2 * Cfg.N, DiffBuffer[1].get());

	// write data
	if (Cfg.datafile)
	{	FILEguard fout(Cfg.datafile, "wt");
		const fftw_real* Up = InBuffer[0].get();
		const fftw_real* Ip = InBuffer[1].get();
		fputs("#t\tU\tI\t∫U\t∫I\tΔU\tΔI\n", fout);
		for (unsigned len = 0; len < Cfg.N; ++len, ++Up, ++Ip)
			fprintf(fout, "%f\t%g\t%g\t%g\t%g\t%g\t%g\n",
				(double)len / Cfg.srate * Cfg.harmonic, *Up, *Ip, Up[Cfg.N], Ip[Cfg.N], Up[2 * Cfg.N], Ip[2 * Cfg.N]);
	}
}

void AnalyzeIn::ExecuteFFT()
{	// forward transformation
	fftwf_execute_r2r(P, InBuffer[0].get(), FFTBuffer[0].get());
	fftwf_execute_r2r(P, InBuffer[1].get(), FFTBuffer[1].get());

	static const double minscale = 1E-15;
	if (Cfg.purgech)
	{	// purge DC (meaningless without DC-coupling)
		FFTBuffer[0].slice(0, Cfg.purgech + 1) *= minscale;
		FFTBuffer[0].slice(Cfg.N - Cfg.purgech, Cfg.purgech) *= minscale;
		FFTBuffer[1].slice(0, Cfg.purgech + 1) *= minscale;
		FFTBuffer[1].slice(Cfg.N - Cfg.purgech, Cfg.purgech) *= minscale;
	}
	// append constant zero to FFT result to simplify analysis
	FFTBuffer[0][Cfg.N] = 0;
	FFTBuffer[1][Cfg.N] = 0;

	// Calculate cross correlation to compensate for constant group delay.
	if (Cfg.crosscorr)
	{	Complex ccv = ExecuteCrossCorrelation(FFTBuffer[0], FFTBuffer[1]);
		// calc linphase to compensate for the group delay
		LinPhase = arg(ccv) / N2f;
	}

	// normalize (including retransformation)
	FFTBuffer[0].slice(0, Cfg.N) *= 1. / Cfg.N;
	FFTBuffer[1].slice(0, Cfg.N) *= 1. / Cfg.N;

	ApplyCalibration();
}

Complex AnalyzeIn::ExecuteCrossCorrelation(const unique_num_array<fftw_real>& in1, const unique_num_array<fftw_real>& in2)
{	assert(in1.size() == in2.size());
	// calculate in1[] * conjugate(in2[])
	// f(0)
	CCBuffer1[0] = in1[0] * in2[0];
	// f(1..N/2-1)
	for (unsigned i = 1; i < Cfg.N / 2; ++i)
	{	Complex c = Complex(in1[i], in1[Cfg.N - i]) * Complex(in2[i], -in2[Cfg.N - i]);
		CCBuffer1[i] = c.real();
		CCBuffer1[Cfg.N - i] = c.imag();
	}
	// f(N/2)
	CCBuffer1[Cfg.N / 2] = in1[Cfg.N / 2] * in2[Cfg.N / 2];

	// do the cross correlation
	fftwf_execute_r2r(PI, CCBuffer1.get(), CCBuffer2.get());

	// use the result
	const double phiinc = M_2PI / Cfg.N;
	Complex r = 0;
	double sum = 0;
	for (unsigned i = 0; i < Cfg.N; ++i)
	{	double amp = CCBuffer2[i];
		amp *= amp;
		//fprintf(stderr, "CC %i\t%g\t%g\n", i, CCBuffer2[i], amp);
		r += polar(amp, phiinc * i);
		sum += amp;
	}
	r /= sum;
	//fprintf(stderr, "CC %g, %g°\n", abs(r), arg(r) * M_180_PI);
	return r;
}
