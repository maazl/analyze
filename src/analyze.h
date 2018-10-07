#ifndef ANALYZE_H_
#define ANALYZE_H_

#include "config.h"
#include "pcmio.h"
#include "mathx.h"
#include "config.h"
#include "utils.h"
#include <fftw3.h>
#include <cmath>
#include <complex>
#include <memory>

#include "unique_array.h"

typedef float fftw_real;


struct filecolumn        ///< configuration reference to a column within a file
{	const char* file;      ///< ... file
	unsigned column;       ///< ... column in file
	constexpr filecolumn() noexcept : file(nullptr), column(1) {}
};

struct action
{	const char* shell;     ///< Command passed to \see system().
	const char* out;       ///< String written to \c stdout.
	constexpr action() noexcept : shell(nullptr), out(nullptr) {}
	void execute() const;
};

struct Config
{	// general options
	bool        help = false;         ///< Show detailed help screen
	const char* cfgout = nullptr;     ///< write effective configuration to this file
	unsigned    srate = 48000;        ///< sampling rate
	bool        floatsamp = false;    ///< read/write floating samples instead of int16
	bool        swapbytes = false;    ///< swap bytes on PCM input
	unsigned    N = 8192;             ///< FFT or analysis block length
	// input options
	const char* infile = nullptr;     ///< PCM input file name
	filecolumn  overwrt[2];           ///< overwrite channel with ...
	const char* rawfile = nullptr;    ///< file name for raw data
	const char* srcfile = nullptr;    ///< file name for input data ASCII output
	double      rref = 1;             ///< value of the reference resistor in impedance measurements
	fftw_real   gainadj[2] = { 1, 1 };///< input gain {l, r}
	unsigned    discsamp = 0;         ///< skip the first samples
	bool        disctrail = false;    ///< consume trailing samples after completion
	bool        diffmode = false;     ///< Differential mode, i.e. denominator I(t) = channel 2 - channel 1
	bool        swapch = false;       ///< Swap input channels L <-> R
	unsigned    addloop = 1;          ///< add raw data before analysis # times
	bool        incremental = false;  ///< incremental mode (add all raw data)
	// output options
	const char* outfile = nullptr;    ///< PCM output file name
	double      outgain = 0.;         ///< output gain
	const char* specfile = nullptr;   ///< file name for frequency domain reference data
	const char* reffile = nullptr;    ///< file name for time domain reference data
	const char* refmode = "w";        ///< open mode for time domain reference data
	double      scalepow = 0;
	// control options
	unsigned    loops = 1;            ///< number of analysis loops
	bool        mpca = false;         ///< analysis method PCA
	bool        mfft = false;         ///< analysis method FFT
	bool        mxy = false;          ///< analysis method XY
	bool        sweep = false;        ///< Use sweep instead of noise
	unsigned    sync = 0;             ///< synchronize cycles before start of measurement
	bool        stereo = false;       ///< Stereo aggregate mode (Toggle harmonics)
	action      init;                 ///< action to take before any processing
	action      plot;                 ///< action to take after analysis step
	action      post;                 ///< action to take after program completion
	// FFT parameter
	const char* datafile = nullptr;   ///< filename for analysis data
	unsigned    winfn = 0;            ///< window function: 0 = rectangle, 1 = Bartlett, 2 = Hanning, 3 = Hamming, 4 = Blackman, 5 = Blackman-Harris
	const char* windowfile = nullptr; ///< file name for window data
	double      fmin = 1E-3;          ///< minimum frequency for FFT analysis
	double      fmax = INFINITY;      ///< minimum frequency for FFT analysis
	double      famin = 1;            ///< ignore frequencies below famin for calculation of screen output
	double      famax = 1E99;         ///< ignore frequencies above famax for calculation of screen output
	double      f_inc = 1;            ///< Absolute increment for harmonic table calculation
	double      f_log = 1;            ///< Relative increment for harmonic table calculation
	unsigned    purgech = 1;          ///< set the first FFT frequencies to 0
	unsigned    binsz = 1;            ///< binsize in FFT channels
	double      fbinsc = 0;           ///< logarithmic binsize: fmin/fmax = 1 + fbinsc
	double      linphase = 0;         ///< linear phase correction [s]
	bool        crosscorr = false;    ///< Calculate and remove time delay by cross correlation
	bool        normalize = false;    ///< normalize L+R to 1. for impedance measurements
	unsigned    harmonic = 0;         ///< analyze up to # harmonics
	double (*weightfn)(double, double, double) = &Config::GetWeight;///< weight function
	// calibration options
	const char* gaininfile = nullptr; ///< file name for gain calibration data
	const char* gainoutfile = nullptr;///< file name for differential gain calibration data
	const char* zeroinfile = nullptr; ///< file name for zero calibration data
	const char* zerooutfile = nullptr;///< file name for differential zero calibration data
	unsigned    lpause = 10;          ///< number of loops between zero calibration parts

	// weight functions
	static double GetWeight(double a1, double a2, double)
	{	double w = 1. / (1. / sqr(a1) + 1. / sqr(a2));
		return std::isfinite(w) ? w : 0.;
	}
	static double GetWeightD(double a1, double a2, double)
	{	double a1q = 1 / sqr(a1);
		double w = 1. / (a1q + sqrt(a1q + 1. / sqr(a1 + a2)) / a2);
		return std::isfinite(w) ? w : 0.;
	}
	static double GetConstWeight(double, double, double)
	{	return 1.;
	}
	static double Get1_fWeight(double, double, double f)
	{	return f ? 1. / f : 0;
	}
};

/// Interface for asynchronous tasks.
struct ITask
{	virtual void Run() = 0;
	virtual ~ITask() {};
};

class AnalyzeIn : public ITask
{private:
	const Config& Cfg;
	const double N2f;                       ///< Frequency bin size

	PCMinput PCMIn;
	FILEguard FIn;
	double LinPhase;
	bool ZeroPart2 = false;                 ///< We are in phase 2 of zero calibration

	// data buffers
	unique_num_array<char> InBufferTmp;     ///< Buffer for raw input
	unique_fftw_arr<fftw_real> InBuffer[2]; ///< Buffer for numerator and denominator input
	unique_fftw_arr<fftw_real> OvrBuffer[2];///< Buffer for overridden numerator/denominator
	unique_fftw_arr<fftw_real> FFTBuffer[2];///< Buffer for FFT(inbuffer[])
	unique_fftw_arr<fftw_real> CCBuffer1;   ///< Buffer for cross correlation temporary data
	unique_fftw_arr<fftw_real> CCBuffer2;   ///< Buffer for cross correlation of outbuffer
	unique_num_array<fftw_real> Window;     ///< Buffer for window function
	unique_num_array<int> Harmonics;        ///< Buffer for harmonics dispatch table

	unique_num_array<double> WSums;         ///< Weights [N+1]
	unique_num_array<Complex> Gain;         ///< Apply gain correction [N / 2 + 1]
	unique_num_array<Complex> GainD;        ///< Result of gain calibration [N / 2 + 1]
	unique_num_array<MatrixC<2,2>> Zero;    ///< Apply matrix correction [N / 2 + 1]
	unique_num_array<MatrixC<2,2>> ZeroD;   ///< Result matrix calibration [N / 2 + 1]

	fftwf_plan P;                           ///< FFT plan for forward transformation
	fftwf_plan PI;                          ///< FFT plan for inverse transformation

 private:
	AnalyzeIn(const Config& cfg);
	/// Setup
	/// @return true if \see Run needs to be called asynchronously.
	bool Setup();
	/// Do the main work
	void Run();
	static void ReadColumn(const filecolumn& src, const unique_fftw_arr<fftw_real>& dst);
	static void CreateWindow(const unique_num_array<fftw_real>& dst, int type);
	void ApplyCalibration(int bin, double f, Complex& U, Complex& I);
	/// Phase unwrapper
	/// Adjusts phase by adding a multiple of 2pi so that the result
	/// is as close as possible to lph.
	/// @param lph last phase
	/// @param  phase current calculated phase
	/// @return unwrapped phase
	static double UnWrap(double lph, double phase)
	{	return !std::isfinite(lph) ? phase : phase - M_2PI * floor((phase - lph) / M_2PI + .5);
	}

	void DoFFT();

 public:
	/// Setup output worker.
	/// @param cfg global configuration.
	/// @return nonzero if there is something to do asynchronously.
	static std::unique_ptr<AnalyzeIn> Setup(const Config& cfg)
	{	std::unique_ptr<AnalyzeIn> ret(new AnalyzeIn(cfg));
		if (!ret->Setup()) ret.reset();
		return ret;
	}

 private:
	class FFTbin
	{public:
		enum StoreRet
		{	BelowMin,    ///< frequency less than fmin
			AboveMax,    ///< frequency above fmax
			Ready,       ///< calculated values available
			Aggregated,  ///< bin used for aggregation only
			Skip         ///< skip this bin because it is a harmonic
		};
	 private:
		struct aggentry
		{	double f;
			double Uabs; ///< Magnitude of numerator
			double Uarg; ///< Phase of numerator
			double Iabs; ///< Magnitude of denominator
			double Iarg; ///< Phase of denominator
			double Zabs; ///< Magnitude of quotient
			double Zarg; ///< Phase of quotient
			double D;    ///< Group delay
			double W;    ///< weight sum
			// internals
			double lf;   ///< last frequency (for numerical derivative)
			double lUarg;///< last phase of numerator
			double lIarg;///< last phase of denominator
			double lZarg;///< last phase of quotient
			double fnext;///< next frequency for bin size
			unsigned binc;///< number of bins accumulated
		};
	 private:
		AnalyzeIn& Parent;

		aggentry agg[2 * HA_MAX + 1];
		int ch;
		aggentry* curagg;
		Complex Zcache;

	 public:
		FFTbin(AnalyzeIn& parent)
		:	Parent(parent)
		{	memset(agg, 0, sizeof agg);
		}
		StoreRet StoreBin(unsigned bin);

		static void PrintHdr(FILE* dst);
		void PrintBin(FILE* dst) const;

		double f() const    { return curagg->f; }    ///< frequency
		Complex U() const   { return std::polar(curagg->Uabs, curagg->Uarg); } ///< voltage, numerator or wanted signal
		double Uabs() const { return curagg->Uabs; } ///< voltage magnitude
		double Uarg() const { return curagg->Uarg; } ///< voltage phase
		Complex I() const   { return std::polar(curagg->Iabs, curagg->Iarg); } ///< current, denominator or reference signal
		double Iabs() const { return curagg->Iabs; } ///< current magnitude
		double Iarg() const { return curagg->Iarg; } ///< current phase
		Complex Z() const   { return Zcache; }       ///< impedance, quotient or relative signal
		double Zabs() const { return curagg->Zabs; } ///< impedance magnitude
		double Zarg() const { return curagg->Zarg; } ///< impedance phase
		double D() const    { return curagg->D / M_2PI; } ///< group delay
		double W() const    { return curagg->W; }    ///< weight
		int    h() const    { return ch; }           ///< harmonic
	};
};

class AnalyzeOut final : public ITask
{private:
	const Config& Cfg;
	const double N2f;     ///< Frequency bin size.
	const unsigned FminI; ///< Smallest used frequency bin.
	const unsigned FmaxI; ///< Highest used frequency bin.
	const double OutLevel;///< Output level (0..1]
	const unsigned LoopCount;///< Number of loops till termination, including zero calibration phases.
 private:
	PCMoutput PCMOut;
	FILEguard FOut;       ///< PCM output file handle or nullptr.
	unsigned FCount;      ///< Number of used frequencies.
	unique_fftw_arr<fftw_real> FFTBuf; ///< Design spectrum
	unique_fftw_arr<char> OutBuf;      ///< PCM data
 private:
	AnalyzeOut(const Config& cfg);
	bool Setup();
	void Normalize(const unique_fftw_arr<fftw_real>& dst);
	virtual void Run();
	/// Play OutBuf LoopCount times or unless termrq.
	void PlayLoop();
 public:
	/// Setup output worker.
	/// @param cfg global configuration.
	/// @return nonzero if there is something to do asynchronously.
	static std::unique_ptr<AnalyzeOut> Setup(const Config& cfg)
	{	std::unique_ptr<AnalyzeOut> ret(new AnalyzeOut(cfg));
		if (!ret->Setup()) ret.reset();
		return ret;
	}
};

#endif // ANALYZE_H_
