#include "analyze.h"
#include "interpolation.h"
#include "mathx.h"
#include "pca.h"
#include "utils.h"

using namespace std;


AnalyzeIn::FFTbin::StoreRet AnalyzeIn::FFTbin::StoreBin(unsigned bin)
{
	curagg = NULL;
	// frequency
	double f = bin * Parent.N2f;
	ch = Parent.Harmonics[bin];
	if ((unsigned)(abs(ch) - 1) >= HA_MAX)
		return Skip;
	curagg = agg + ch + HA_MAX;
	unsigned base = bin;        //ch != 0 ? bin/abs(ch) : bin;
	// retrieve coefficients
	Complex U(Parent.FFTBuffer[0][bin], bin && bin != Parent.Cfg.N / 2 ? Parent.FFTBuffer[0][Parent.Cfg.N - bin] : 0);
	Complex I(Parent.FFTBuffer[1][base], bin && bin != Parent.Cfg.N / 2 ? Parent.FFTBuffer[1][Parent.Cfg.N - base] : 0);
	// calibration
	Parent.ApplyCalibration(bin, f, U, I);
	// phase correction
	U *= Complex(cos(Parent.LinPhase * f), sin(Parent.LinPhase * f));
	// calc Y
	Complex Z(U / I);
	// convert to polar
	double Uabs = abs(U);
	double Uarg = UnWrap(curagg->lUarg, arg(U));
	double Iabs = abs(I);
	double Iarg = UnWrap(curagg->lIarg, arg(I));
	double Zabs = abs(Z);
	double Zarg = UnWrap(curagg->lZarg, arg(Z));
	// group delay
	double D = (Zarg - curagg->lZarg) / (f - curagg->lf);
	// store values for next point
	curagg->lUarg = Uarg;
	curagg->lIarg = Iarg;
	curagg->lZarg = Zarg;
	curagg->lf = f;

	// weight
	double w = Parent.Cfg.weightfn(Uabs, Iabs, f);
	if (curagg->binc == 0)
	{	// init
		curagg->f = f * w;
		curagg->Uabs = Uabs * w;
		curagg->Uarg = Uarg * w;
		curagg->Iabs = Iabs * w;
		curagg->Iarg = Iarg * w;
		curagg->Zabs = Zabs * w;
		curagg->Zarg = Zarg * w;
		curagg->D = D * w;
		curagg->W = w;
		curagg->fnext = f * (1 + Parent.Cfg.fbinsc) - Parent.N2f;
	} else
	{
		curagg->f += f * w;
		curagg->Uabs += Uabs * w;
		curagg->Uarg += Uarg * w;
		curagg->Iabs += Iabs * w;
		curagg->Iarg += Iarg * w;
		curagg->Zabs += Zabs * w;
		curagg->Zarg += Zarg * w;
		curagg->D += D * w;
		curagg->W += w;
	}
	++curagg->binc;
	if (f < curagg->fnext)
		return Aggregated;
	w = curagg->W;
	curagg->f /= w;
	curagg->Uabs /= w;
	curagg->Uarg /= w;
	curagg->Iabs /= w;
	curagg->Iarg /= w;
	curagg->Zabs /= w;
	curagg->Zarg /= w;
	curagg->D /= w;
	curagg->binc = 0;
	Zcache = polar(curagg->Zabs, curagg->Zarg);
	/*if (curagg->f < fmin)
	 return BelowMin;
	 if (curagg->f > fmax)
	 return AboveMax;*/
	return Ready;
}

void AnalyzeIn::FFTbin::PrintHdr(FILE* dst)
{
	fputs("#f\t|U|\targ U\t|I|\targ I\t|Z|\targ Z\tZ real\tZ imag\tweight\tdelay\tharmon.\n", dst);
}

void AnalyzeIn::FFTbin::PrintBin(FILE* dst) const
{
	fprintf(dst, "%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%6i\n",
	// f    |Hl|      phil               |Hr|      phir
	   f(), Uabs(), Uarg() * M_180_PI, Iabs(), Iarg() * M_180_PI,
	// |Hl|/|Hr| phil-phir          re          im
	   Zabs(), Zarg() * M_180_PI, Z().real(), Z().imag(),
	// weight delay harmonic
	   W(), D(), h());
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

void AnalyzeIn::ApplyCalibration(int bin, double f, Complex& U, Complex& I)
{
	if (Cfg.gaininfile && !Cfg.gainoutfile)
		U /= Gain[bin];
	if (Cfg.gainoutfile)
	{	double weight = sqrt(Cfg.weightfn(abs(U), abs(I), f));
		WSums[bin] += weight;
		U /= Gain[bin];
		GainD[bin] += I / U * weight;
	}
	if (Cfg.normalize)
	{	Complex s = 1. / (U + I);
		U *= s;
		I *= s;
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
}

void AnalyzeIn::DoFFT()
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
	{	// calculate outbuffer1[] * conjugate(outbuffer2[])
		// f(0)
		CCBuffer1[0] = FFTBuffer[0][0] * FFTBuffer[1][0];
		// f(1..N/2-1)
		for (unsigned i = 1; i < Cfg.N / 2; ++i)
		{	Complex c = Complex(FFTBuffer[0][i], FFTBuffer[0][Cfg.N - i]) * Complex(FFTBuffer[1][i], FFTBuffer[1][Cfg.N - i]);
			CCBuffer1[i] = c.real();
			CCBuffer1[Cfg.N - i] = c.imag();
		}
		// f(N/2)
		CCBuffer1[Cfg.N / 2] = FFTBuffer[0][Cfg.N / 2] * FFTBuffer[1][Cfg.N / 2];

		// do the cross correlation
		fftwf_execute_r2r(PI, CCBuffer1.get(), CCBuffer2.get());

		// use the result
		double phiinc = M_2PI / Cfg.N;
		double asum = 0;
		double bsum = 0;
		double sum = 0;
		for (unsigned i = 0; i < Cfg.N; ++i)
		{
			double amp = CCBuffer2[i];
			amp *= amp;
			asum += cos(phiinc * i) * amp;
			bsum += sin(phiinc * i) * amp;
			sum += amp;
		}
		asum /= sum;
		bsum /= sum;

		// calc linphase to compensate for the group delay
		LinPhase = atan2(bsum, asum) / N2f;
	}
}

AnalyzeIn::AnalyzeIn(const Config& cfg)
	// setup input
:	Cfg(cfg)
,	N2f((double)Cfg.srate / Cfg.N)
,	PCMIn(Cfg.floatsamp ? Format::F32 : Cfg.swapbytes ? Format::I16_SWAP : Format::I16, Cfg.diffmode, &Cfg.gainadj)
,	LinPhase(Cfg.linphase * M_2PI)
{	// allocate buffers
	InBufferTmp.reset(PCMIn.BytesPerSample * Cfg.N);
	Window.reset(Cfg.N);
	if (Cfg.mxy)
	{	// reserve space for integrals and differentials too
		InBuffer[0].reset(3 * Cfg.N);
		InBuffer[1].reset(3 * Cfg.N);
		FFTBuffer[0].reset(3 * Cfg.N + 1);
		FFTBuffer[1].reset(3 * Cfg.N + 1);
	} else
	{	InBuffer[0].reset(Cfg.N);
		InBuffer[1].reset(Cfg.N);
		FFTBuffer[0].reset(Cfg.N + 1);
		FFTBuffer[1].reset(Cfg.N + 1);
	}
	if (Cfg.crosscorr && Cfg.mpca)
	{	CCBuffer1.reset(Cfg.N);
		CCBuffer2.reset(Cfg.N);
	}
	WSums.reset(Cfg.N / 2 + 1);
	WSums.clear();
	Harmonics.reset(Cfg.N / 2 + 1);
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
	// create harmonics table
	{	Harmonics.clear();
		int sign = 1;
		for (unsigned i = (unsigned)ceil(Cfg.fmin / N2f); i <= floor(fmax / N2f); ++i)
		{	if (i)
			{	for (unsigned j = 1; j <= Cfg.harmonic && i * j <= Cfg.N / 2; ++j)
					if (Harmonics[i * j])
						goto next_f;
				for (unsigned j = 1; i * j <= Cfg.N / 2; ++j)
					Harmonics[i * j] = j * sign;
			}
			if (Cfg.stereo)
				sign = -sign;
			i = (unsigned)floor(i * Cfg.f_log + Cfg.f_inc - .5);
		 next_f:;
		}
		/*FILE* f = fopen("harm.dat", "w");
		 for (unsigned i = 0; i <= N/2; ++i)
		 fprintf(f, "%12g %8i\n", i*freq/N, harmonics[i]);
		 fclose(f);*/
	}

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
	{	PolarFileInterpolation ip(Cfg.gaininfile, 1);
		unsigned i = 0;
		for (auto& g : Gain)
		{	auto& row = ip.Get(i++ * N2f);
			g = Complex(row[1], row[2]);
		}
	}
	// read zero file
	if (Cfg.zeroinfile)
	{	PolarFileInterpolation ip(Cfg.zeroinfile, 4);
		unsigned i = 0;
		for (auto& z : Zero)
		{	auto& row = ip.Get(i++ * N2f);
			z[0][0] = Complex(row[1], row[2]);
			z[0][1] = Complex(row[3], row[4]);
			z[1][0] = Complex(row[5], row[6]);
			z[1][1] = Complex(row[7], row[8]);
			z = inverse(z); // invert matrix to save time during analysis
		}
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

	const unique_num_array<fftw_real> input[2] =
	{	unique_num_array<fftw_real>(InBuffer[0].slice(0, Cfg.N))
	,	unique_num_array<fftw_real>(InBuffer[1].slice(0, Cfg.N))
	};
	input[0].clear(); // init with 0 because of incremental mode
	input[1].clear();

 restart_zero:
	// operation loop
	unsigned loop = Cfg.loops;
	unsigned addloops = 0;
	do
	{
		if (FIn)
		{	// read data
			fread2(InBufferTmp.get(), InBufferTmp.size(), FIn);
			// write raw data
			if (Cfg.rawfile)
				PCMIn.ASCIIdump(FILEguard(Cfg.rawfile, "wt"), InBufferTmp);

			// reset min/max
			PCMIn.reset();
			// convert data
			PCMIn.convert(input[Cfg.swapch], input[!Cfg.swapch], InBufferTmp, Cfg.incremental || addloops);
			// write raw status
			fprintf(stderr, "\nmin:\t%f\t%f\nmax:\t%f\t%f\n",
				PCMIn.Limits[0].Min, PCMIn.Limits[1].Min, PCMIn.Limits[0].Max, PCMIn.Limits[1].Max);
		}

		if (Cfg.incremental || addloops)
		{	if (OvrBuffer[0])
				input[0] += OvrBuffer[0];
			if (OvrBuffer[1])
				input[1] += OvrBuffer[1];
		} else
		{	if (OvrBuffer[0])
				input[0] = OvrBuffer[0];
			if (OvrBuffer[1])
				input[1] = OvrBuffer[1];
		}

		if (++addloops < Cfg.addloop)
			continue; // add more data

		// write source data
		if (Cfg.srcfile)
		{	FILEguard out(Cfg.srcfile, "wt");
			size_t len = InBuffer[0].size();
			const fftw_real* sp1 = InBuffer[0].get();
			const fftw_real* sp2 = InBuffer[1].get();
			while (len--)
				fprintf(out, "%g\t%g\n", *sp1++, *sp2++);
		}

		if (Cfg.winfn && (Cfg.incremental || addloops))
		{	// optimization: apply window function later
			input[0] *= Window;
			input[1] *= Window;
		}
		addloops = 0;

		if (Cfg.mfft & Cfg.mpca)
		{	// FFT
			DoFFT();

			FFTBuffer[0] *= sqrt(1. / Cfg.N);
			FFTBuffer[1] *= sqrt(1. / Cfg.N);

			// write data
			FILEguard tout = NULL;
			if (Cfg.datafile)
			{	tout = checkedopen(Cfg.datafile, "wt");
				FFTbin::PrintHdr(tout);
			}

			PCA<2> pcaRe;
			PCA<3> pcaIm;
			double PCAdataRe[2];
			double PCAdataIm[3];
			// some values are const
			PCAdataRe[1] = 1;
			//PCAdataIm[1] = 1;

			FFTbin calc(*this);

			// 1st line
			for (size_t len = 0; len <= Cfg.N / 2; ++len)
			{	// do calculations and aggregations
				switch (calc.StoreBin(len))
				{case FFTbin::AboveMax:
					// write
					if (tout && calc.h() > 1)
						calc.PrintBin(tout);
				 default:
					//case FFTbin::BelowMin:
					//case FFTbin::Aggregated:
					//case FFTbin::Skip:
					continue;
				 case FFTbin::Ready:
					// write
					if (tout)
						calc.PrintBin(tout);
				}

				if (calc.f() < Cfg.famin || calc.f() > Cfg.famax)
					continue;
				// component analysis
				PCAdataRe[0] = calc.Z().real();
				//PCAdataRe[2] = 1/af;
				//PCAdataRe[3] = f;
				PCAdataIm[0] = calc.Z().imag(); // fit imaginary part in conductivity
				PCAdataIm[1] = 1 / calc.f();
				PCAdataIm[2] = calc.f();
				//PCAdataIm[3] = 1/af;
				//fprintf(stderr, "Re: %12g %12g %12g %12g\n", calc.f(), PCAdataRe[0], PCAdataRe[1], calc.W());

				// add values
				pcaRe.Store(PCAdataRe, calc.W());
				pcaIm.Store(PCAdataIm, calc.W());
			}

			// calculate summary
			VectorD<1> resRe = pcaRe.Result();
			VectorD<2> resIm = pcaIm.Result();

			fprintf(stderr, "resRe: %12g %12g %12g\n", resRe[0], resRe[1], resRe[2]);
			//printf("resRe: %12g %12g %12g %12g\n", resRe[0], resRe[1], resRe[2], resRe[3]);
			fprintf(stderr, "resIm: %12g %12g %12g %12g\n", resIm[0], resIm[1], resIm[2], resIm[3]);

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
		else if (Cfg.mpca)
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
		else if (Cfg.mfft)
		{
			DoFFT();

			FFTBuffer[0].slice(0, Cfg.N) *= sqrt(1. / Cfg.N);
			FFTBuffer[1].slice(0, Cfg.N) *= sqrt(1. / Cfg.N);

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
				FFTbin::PrintHdr(tout);
			}

			FFTbin calc(*this);

			for (size_t len = 0; len <= Cfg.N / 2; ++len)
			{	// do calculations and aggregations
				switch (calc.StoreBin(len))
				{case FFTbin::AboveMax:
					// write
					if (tout && abs(calc.h()) > 1)
						calc.PrintBin(tout);
				 default:
					//case FFTbin::BelowMin:
					//case FFTbin::Aggregated:
					//case FFTbin::Skip:
					continue;
				 case FFTbin::Ready:
					// write
					if (tout)
						calc.PrintBin(tout);
				}

				if (calc.f() < Cfg.famin || calc.f() > Cfg.famax)
					continue;
				//fprintf(stderr, "FW %12g %12g\n", calc.f(), calc.W());
				// average
				++nsum;
				wsum += calc.W();
				// resistivity
				Rsum += calc.W() * calc.Z().real();
				R2sum += calc.W() * sqr(calc.Z().real());
				// L & C
				L2sum += calc.W() * sqr(calc.f());
				// LCsum += weight; == wsum
				C2sum += calc.W() / sqr(calc.f());
				Lsum += calc.W() * calc.Z().imag() * calc.f();
				Csum += calc.W() * calc.Z().imag() / calc.f();
				d2sum = calc.W() * sqr(calc.Z().imag());
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
				"\nreal: R [Ω]\t%12g ± %g\n"
				"imag.: C [µF]\t%12g ± %g\n"
				"imag.: L [µH]\t%12g ± %g\n",
				Cfg.rref * R, Cfg.rref * RE, 1E6 / (Cfg.rref * C * M_2PI), 1E6 / (CE * Cfg.rref * M_2PI),
				1E6 * Cfg.rref * L / M_2PI, 1E6 * Cfg.rref * LE / M_2PI);
			if (Cfg.crosscorr)
				fprintf(stderr, "delay [ms]\t%12g\n", 1E3 / M_2PI * LinPhase);
		}
		else if (Cfg.mxy)
		{	// we need to do an FFT here, at least for the calibration
			DoFFT();
			// normalize (including retransformation)
			FFTBuffer[0].slice(0, Cfg.N) *= sqrt(1. / sqr(Cfg.N));
			FFTBuffer[1].slice(0, Cfg.N) *= sqrt(1. / sqr(Cfg.N));

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
				// calibration
				ApplyCalibration(len, f, U, I);
				// store data
				*a1 = U.real();
				*b1 = U.imag();
				*a2 = I.real();
				*b2 = I.imag();
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
			fftwf_execute_r2r(PI, FFTBuffer[0].get() + Cfg.N, InBuffer[0].get() + Cfg.N);
			fftwf_execute_r2r(PI, FFTBuffer[0].get() + 2 * Cfg.N, InBuffer[0].get() + 2 * Cfg.N);
			fftwf_execute_r2r(PI, FFTBuffer[1].get(), InBuffer[1].get());
			fftwf_execute_r2r(PI, FFTBuffer[1].get() + Cfg.N, InBuffer[1].get() + Cfg.N);
			fftwf_execute_r2r(PI, FFTBuffer[1].get() + 2 * Cfg.N, InBuffer[1].get() + 2 * Cfg.N);

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

		Cfg.plot.execute();

		// undo window function because of incremental mode
		// TODO: this causes a loss of precision!
		if (Cfg.winfn && Cfg.incremental)
		{	input[0] /= Window;
			input[1] /= Window;
		}

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
			for (int loop = Cfg.lpause; loop;)
			{	fprintf(stderr, "\r%u ", loop);
				PCMIn.discard(FIn, Cfg.N, InBufferTmp);
				--loop;
			}
			puts("\nNow at part 2.");
			goto restart_zero;
		} else // part 2
		{	FILEguard fz(Cfg.zerooutfile, "wb");
			fputs("#f\tU->U re\tU->U im\tU->I re\tU->I im\tI->U re\tI->U im\tI->I re\tI->I im\t"
				"abs U->U\targ U->U\tabs U->I\targ U->I\tabs I->U\targ I->U\tabs I->I\targ I->I\n", fz);
			unsigned i = 0;
			for (const auto& z : ZeroD)
			{	// scale to fit det z == 1
				Complex sca = 1. / sqrt(det(z));
				if (((z[0][0] + z[1][1]) * sca).real() < 0)
					sca = -sca;
				auto z_ = z * sca;
				fprintf(fz, "%f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					i++ * N2f,
					z_[0][0].real(), z_[0][0].imag(), z_[0][1].real(), z_[0][1].imag(),
					z_[1][0].real(), z_[1][0].imag(), z_[1][1].real(), z_[1][1].imag(),
					abs(z_[0][0]), arg(z_[0][0]) * M_180_PI, abs(z_[0][1]), arg(z_[0][1]) * M_180_PI,
					abs(z_[1][0]), arg(z_[1][0]) * M_180_PI, abs(z_[1][1]), arg(z_[1][1]) * M_180_PI);
				/*fprintf(fz, "%f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					i++ * N2f,
					z[0][0].real(), z[0][0].imag(), z[0][1].real(), z[0][1].imag(),
					z[1][0].real(), z[1][0].imag(), z[1][1].real(), z[1][1].imag(),
					abs(z[0][0]), arg(z[0][0]) * M_180_PI, abs(z[0][1]), arg(z[0][1]) * M_180_PI,
					abs(z[1][0]), arg(z[1][0]) * M_180_PI, abs(z[1][1]), arg(z[1][1]) * M_180_PI);*/
			}
		}
	}

	fputs("completed.", stderr);

	Cfg.post.execute();

	// read until EOF
	if (Cfg.disctrail && !termrq)
		while (fread(InBufferTmp.get(), sizeof InBufferTmp[0], InBufferTmp.size(), FIn) > 0);
}
