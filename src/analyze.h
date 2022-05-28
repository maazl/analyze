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
#include "filereader.h"

typedef float fftw_real;

typedef double (*WeightFn)(double, double, double);

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
	double      srate = 48000.;       ///< sampling rate
	Format      format = Format::I16;///< read/write samples format
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
	bool        minmax = false;       ///< Show min max
	unsigned    addloop = 1;          ///< add raw data before analysis # times
	bool        incremental = false;  ///< incremental mode (add all raw data)
	// output options
	const char* outfile = nullptr;    ///< PCM output file name
	double      outgain = 0.;         ///< output gain
	bool        nohdr = false;        ///< no WAV header
	unsigned    outch = 0;            ///< output channel
	const char* rspecfile = nullptr;  ///< file name for frequency domain reference data input
	const char* specfile = nullptr;   ///< file name for frequency domain reference data output
	const char* reffile = nullptr;    ///< file name for time domain reference data
	const char* refmode = "w";        ///< open mode for time domain reference data
	double      scalepow = 0;         ///< exponent of f for power distribution
	double      chirp = 0;            ///< Use chirp phase with given duty cycle; 0 = random phase
	double      smooth = 0;           ///< Smoothen the range of used frequencies [in FFT channels]
	// control options
	unsigned    loops = 1;            ///< number of analysis loops, 0 = infinite
	bool        mpca = false;         ///< analysis method PCA
	bool        mfft = false;         ///< analysis method FFT
	bool        mxy = false;          ///< analysis method XY
	bool        sweep = false;        ///< Use sweep instead of noise
	bool        stereo = false;       ///< Stereo aggregate mode (Toggle harmonics)
	action      setup;                ///< action to take at program startup
	action      init;                 ///< action to take before any data processing
	action      plot;                 ///< action to take after analysis step
	action      post;                 ///< action to take after program completion
	// FFT parameter
	const char* datafile = nullptr;   ///< filename for analysis data
	const char* irfile = nullptr;     ///< filename for impulse response data
	unsigned    winfn = 0;            ///< window function: 0 = rectangle, 1 = Bartlett, 2 = Hanning, 3 = Hamming, 4 = Blackman, 5 = Blackman-Harris
	const char* windowfile = nullptr; ///< file name for window data
	double      fmin = 1E-3;          ///< minimum frequency for FFT analysis
	double      fmax = INFINITY;      ///< minimum frequency for FFT analysis
	double      famin = 1;            ///< ignore frequencies below famin for calculation of screen output
	double      famax = 1E99;         ///< ignore frequencies above famax for calculation of screen output
	double      f_inc = 1;            ///< Absolute increment for harmonic table calculation
	double      f_log = 1;            ///< Relative increment for harmonic table calculation
	unsigned    purgech = 1;          ///< set the first FFT frequencies to 0
	double      binsz = 1;            ///< binsize in FFT channels
	double      fbinsc = 1;           ///< logarithmic binsize: fmin/fmax = 1 + fbinsc
	double      linphase = 0;         ///< linear phase correction [s]
	bool        crosscorr = false;    ///< Calculate and remove time delay by cross correlation
	unsigned    harmonic = 1;         ///< analyze up to # harmonics
	WeightFn    weightfn = &Config::GetWeight;///< weight function
	// synchronization options
	unsigned    sync = 0;             ///< synchronize cycles before start of measurement
	unsigned    syncch = 3;           ///< sync channel 1 = left, 2 = right, 3 = both
	//bool        syncfall = true;      ///< use all frequencies in the interval [famin, famax] rather than the design function
	double      synclevel = .2;       ///< minimum SNR ratio of successful synchronization compared to the SNR of the autocorrelation function
	double      syncend = .8;         ///< decrease of cross correlation level that is identified as end of sync.
	double      predelay = .95;       ///< gap before a measurement starts or after a frequency change in units of FFT cycles (N).
	// calibration options
	const char* gaininfile = nullptr; ///< file name for gain calibration data
	const char* gainoutfile = nullptr;///< file name for differential gain calibration data
	const char* zeroinfile = nullptr; ///< file name for zero calibration data
	const char* zerooutfile = nullptr;///< file name for differential zero calibration data
	bool        normalize = false;    ///< normalize L+R to 1. for two point calibration
	double      lpause = 10;          ///< number of seconds between zero calibration parts
	unsigned    PauseLoops() const    ///< Number of loops for pause between zero calibration parts.
	{	return (unsigned)ceil(lpause / N * srate); }

	// weight functions
	static double GetWeight1(double a1, double a2, double)
	{	return a1;
	}
	static double GetWeight2(double a1, double a2, double)
	{	return a2;
	}
	static double GetWeight(double a1, double a2, double)
	{	double w = 1. / (1. / sqr(a1) + 1. / sqr(a2));
		if (!std::isfinite(w))
			w = 0;
		return w + 1E-99;
	}
	static double GetWeightD(double a1, double a2, double)
	{	double a1q = 1 / sqr(a1);
		double w = 1. / (a1q + sqrt(a1q + 1. / sqr(a1 + a2)) / a2);
		if (!std::isfinite(w))
			w = 0;
		return w + 1E-99;
	}
	static double GetConstWeight(double, double, double)
	{	return 1.;
	}
	static double Get1_fWeight(double, double, double f)
	{	return f ? 1. / f : 1E-99;
	}
};

/// Safe wrapper around \c fftwf_plan, calls \c fftwf_destroy_plan on destruction.
class unique_fftwf_plan final
{	fftwf_plan Plan;
 public:
	unique_fftwf_plan() : Plan(NULL) {}
	unique_fftwf_plan(fftwf_plan plan) : Plan(plan) {}
	unique_fftwf_plan(unique_fftwf_plan&) = delete;
	~unique_fftwf_plan() { fftwf_destroy_plan(Plan); }
	operator fftwf_plan() const { return Plan; }
	unique_fftwf_plan& operator=(fftwf_plan plan) { fftwf_destroy_plan(Plan); Plan = plan; return *this; }
	unique_fftwf_plan& operator=(unique_fftwf_plan&) = delete;
	fftwf_plan get() const { return Plan; }
};

/// Interface for asynchronous tasks.
struct ITask
{	virtual void Run() = 0;
	virtual ~ITask() {};
};

struct SetupData
{
	/// Design spectrum in packed half complex format, i.e [a(0), a(1) ... a(N/2), b(N/2-1) ... b(1)].
	/// b(N/2) and b(0) are zero by definition of half complex FFT.
	unique_fftw_arr<fftw_real> Design;
	/// Used frequencies and harmonics.
	/// @details The table is of size N/2+1. Each entry corresponds to the matching FFT frequency bin
	/// and can have one of the following values:
	/// - 0 := the frequency is unused.
	/// - 1 := the frequency is directly used.
	/// - n := the frequency is the n-th harmonic of a directly used frequency and not exported.
	/// - < 0 := the frequency is used by the second channel rather than the first one.
	unique_num_array<int> Harmonics;
	/// Number of used frequencies.
	unsigned FCount;
};

/// Analysis part
class AnalyzeIn : public ITask
{private:
	const Config& Cfg;                      ///< global configuration
	const SetupData& SD;                    ///< setup data from AnalyzeOut
	const double N2f;                       ///< Frequency bin size, i.e. sampling rate / N

	PCMinput PCMIn;
	FILEguard FIn;                          ///< Input data stream
	double LinPhase;                        ///< Delay between signal and reference in seconds/2Pi
	double SNRAC;                           ///< Signal to noise of autocorrelation function
	unsigned char ZeroPart = 0;             ///< Part of matrix calibration, count down
	long FileStart = 0;

	bool NeedSync;                          ///< True if Synchronization is (still) required before main data processing.
	unsigned NextSamples;                   ///< Number of samples to read for the next loop. Typically Cfg.N.
	Complex SyncPhaseVector;                ///< Aggregated result from ExecuteCrossCorrelation
	double SyncWeight;                      ///< Sync cycle counter, < 0 => not synced since # cycles, > 0 sync since # cycles.

	// data buffers
	unique_num_array<char> InBufferTmp;     ///< Buffer for raw input
	unique_fftw_arr<fftw_real> InBuffer[2]; ///< Buffer for numerator and denominator input
	unique_fftw_arr<fftw_real> IntBuffer[2];///< Buffer for integral (XY mode only)
	unique_fftw_arr<fftw_real> DiffBuffer[2];///< Buffer for differential (XY mode only)
	unique_fftw_arr<fftw_real> OvrBuffer[2];///< Buffer for overridden numerator/denominator
	unique_fftw_arr<fftw_real> FFTBuffer[2];///< Buffer for FFT(inbuffer[])
	unique_fftw_arr<fftw_real> TmpBuffer[3];///< Buffers for temporary data
	unique_num_array<fftw_real> Window;     ///< Buffer for window function

	unique_num_array<fftw_real> Input[2];   ///< First Cfg.N samples slice of InBuffer
	unique_num_array<double> WSums;         ///< Weights [N+1]
	unique_num_array<Complex> Gain;         ///< Apply gain correction [N / 2 + 1]
	unique_num_array<Complex> GainD;        ///< Result of gain calibration [N / 2 + 1]
	unique_num_array<MatrixC<2,2>> Zero;    ///< Apply matrix correction [N / 2 + 1]
	unique_num_array<MatrixC<3,2>> ZeroD;   ///< Result matrix calibration [N / 2 + 1], [ZeroPart]

	unique_fftwf_plan P;                    ///< FFT plan for forward transformation
	unique_fftwf_plan PI;                   ///< FFT plan for inverse transformation

	class FFTWorker;
	class ImpulseResponseWorker;

	std::unique_ptr<ImpulseResponseWorker> IRWorker[2];

 private:
	/// Do the main work
	void Run();
	/// Synchronize with sync pattern
	/// @post NextSamples (out) [1..Cfg.N] number of samples requested for the next sync try.
	/// If the result is less than Cfg.N the remaining samples should stay in the input buffer.
	/// @p NeedSync true: further synchronization needed; false: synchronization succeeded, start main measurement.
	void DoSync();
	/// Analyze samples in InBuffer by PCA analysis
	void DoPCA();
	/// Analyze samples in InBuffer by FFT analysis
	void DoFFT(bool ch2);
	/// Analyze samples in InBuffer by FFT & PCA analysis
	void DoFFTPCA(bool ch2);
	/// Analyze samples in InBuffer by hysteresis analysis
	void DoXY();

	void FinishImpulseResponse();

	void ExecuteFFT();
	/// Calculate cross correlation
	/// @param in1 FFT of the first input.
	/// @param in2 FFT of the second input.
	/// @return Relative phase vector between first and second input in units of the base frequency.
	/// The amplitude corresponds to the SNR of the cross correlation. The higher the less noise.
	/// But keep in mind that the autocorrelation of the reference signal also adds noise.
	Complex ExecuteCrossCorrelation(const unique_num_array<fftw_real>& in1, const unique_num_array<fftw_real>& in2);

	void PrintHdr(FILE* dst) const;
	void PrintBin(FILE* dst, const FFTWorker& fft) const;

	static void ReadColumn(const filecolumn& src, const unique_fftw_arr<fftw_real>& dst);
	static void CreateWindow(const unique_num_array<fftw_real>& dst, int type);
	void ApplyCalibration();

 public:
	/// Create input worker.
	/// @param cfg global configuration.
	/// @param sd setup information, created by the output worker.
	AnalyzeIn(const Config& cfg, const SetupData& sd);
	/// Setup input worker.
	/// @return true if \see Run needs to be called asynchronously.
	bool Setup();

 private:
	/// Phase unwrapper
	/// Adjusts phase by adding a multiple of 2π
	/// so that the result is as close as possible to the last result.
	struct Unwrapper
	{	double Phase;  ///< Unwrapped phase
		Unwrapper() : Phase(0) {}
		/// Unwrap phase
		/// @param phase current calculated phase
		/// @return unwrapped phase
		double Unwrap(double phase)
		{	return Phase = phase - M_2PI * floor((phase - Phase) / M_2PI + .5);
		}
	};

	struct GroupDelay : public Unwrapper
	{	double Delay;
		/// Unwrap phase and calculate gropup delay.
		/// The unwrapped phase can be extracted by operator double() after the call.
		/// @param phase current calculated phase
		/// @param deltaf frequency difference to the last phase
		/// @return group delay
		double Unwrap(double phase, double deltaf)
		{	double lph = Phase;
			/*Unwrapper::Unwrap(phase);
			fprintf(stderr, "lph = %f, ph = %f, df = %f => ph = %f, delay = %f\n", lph, phase, deltaf, Phase, (Phase - lph) / deltaf);
			return Delay = (Phase - lph) / deltaf;*/
			return Delay = (Unwrapper::Unwrap(phase) - lph) / deltaf;
		}
	};

	struct AggEntry
	{	struct VE
		{	double Uabs; ///< Magnitude of numerator
			double Uarg; ///< Phase of numerator
			double Zabs; ///< Magnitude of quotient
			double Zarg; ///< Phase of quotient
			double D;    ///< Group delay
			double W;    ///< weight sum
		};
		double F;
		double Iabs;   ///< Magnitude of denominator
		double Iarg;   ///< Phase of denominator
		VE Val[HA_MAX];///< Values per harmonic
		AggEntry() { memset(this, 0, sizeof *this); }
		void next() { memset(this, 0, (char*)&lF - (char*)this); }
		void add(double f, Complex I, Complex* U, unsigned harm, WeightFn wfn);
		void finish(unsigned harm);
	 private:
		double lF;     ///< last frequency (for numerical derivative)
		Unwrapper lIarg;///< last phase of denominator
		Unwrapper lUarg[HA_MAX];///< last phase of numerator
		GroupDelay lZarg[HA_MAX];///< last phase of quotient
	};

	class FFTWorker
	{public:
		enum StoreRet
		{	BelowMin,    ///< frequency less than fmin
			AboveMax,    ///< frequency above fmax
			Aggregated,  ///< bin used for aggregation only
			Ready        ///< calculated values available
		};
	 private:
		struct FFTAgg : AggEntry
		{	unsigned Bins;///< Number of bins accumulated
		 	unsigned Harm;///< number of used entries in Val
			double NextF;///< Next frequency
			FFTAgg() : Bins(0), Harm(0), NextF(0) {}
		};

		AnalyzeIn& Parent;

		int Ch;        ///< Channel, i.e. index to Agg
		FFTAgg Agg[2]; ///< Aggregation per channel

	 public:
		FFTWorker(AnalyzeIn& parent) : Parent(parent) {}
		/// Aggregate FFT bin
		/// @param bin FFT bin, [0, Cfg.N/2]
		/// @param ch Channel, 0 for 1st channel.
		StoreRet StoreBin(unsigned bin, int ch);
		/// Current channel
		int ch() const { return Ch; }
		/// Retrieve aggregated information for harmonic of current frequency and channel
		const AggEntry& ret() const { return Agg[Ch]; }
	};

	class SweepWorker
	{private:
		struct SweepBin
		{	Complex Xi;
			double Sum;
			double Sum2;
			SweepBin() : Sum(0), Sum2(0) {}
			void store(Complex vec, double val)
			{	Xi += vec * val;
				Sum += val;
				Sum2 += sqr(val);
			}
			void finish(unsigned n)
			{	Xi /= n * M_SQRT1_2;
				Sum /= n;
				Sum2 /= n;
			}
			double noise() const { return sqrt(Sum2 - norm(Xi) - sqr(Sum)); }
		};

	 private:
		AnalyzeIn& Parent;
		FILEguard Out;           ///< Target data file
		unsigned SweepFrequency; ///< Current frequency index
		unsigned Channel;        ///< Current channel
		SweepBin Ana[1 + HA_MAX];///< Reference & per harmonic
		GroupDelay Delay[2];     ///< per channel group delay
		Unwrapper Phase[2][HA_MAX-1];///< per channel & harmonic

		unsigned LastFrequency;

	 public:
		SweepWorker(AnalyzeIn& parent);
		void DoSweep(unsigned block);
	 private:
		SweepWorker(const SweepWorker&) = delete;
		void operator=(const SweepWorker&) = delete;
		void Print();
	};

	class ImpulseResponseWorker
	{private:
		const AnalyzeIn& Parent;
		const unique_num_array<fftw_real>& FD;
		const unique_num_array<fftw_real>& TD;
		PolarInterpolation Interpolation;
		unsigned NextBin;
	 private:
		void Process();
		ImpulseResponseWorker(const ImpulseResponseWorker&) = delete;
		void operator=(const ImpulseResponseWorker&) = delete;
	 public:
		/// Create impulse response worker
		/// @param fd Place frequency domain data here
		/// @param td Place time domain result here
		ImpulseResponseWorker(AnalyzeIn& parent, const unique_num_array<fftw_real>& fd, const unique_num_array<fftw_real>& td);
		void Reset();
		void Feed(double f, Complex value);
		void Finish();
	};
};

/// Reference signal part
class AnalyzeOut final : public ITask, public SetupData
{private:
	const Config& Cfg;    ///< global configuration
	const double N2f;     ///< Frequency bin size.
	const double MinQuant;///< Energy quantile function at fmin.
	const double MaxQuant;///< Energy quantile function at fmax.
	const double FminSmo; ///< fmin where smooting ends.
	const double FmaxSmo; ///< fmax where smooting starts.
	const unsigned FminI; ///< Smallest used frequency bin.
	const unsigned FmaxI; ///< Highest used frequency bin.
	const double OutLevel;///< Output level (0..1]
	const unsigned LoopCount;///< Number of loops till termination, including zero calibration phases.
 private:
	PCMoutput PCMOut;
	FILEguard FOut;       ///< PCM output file handle or nullptr.
	unique_fftw_arr<char> OutBuf;      ///< PCM data
 private:
	static fftw_real MaxAbs(const unique_fftw_arr<fftw_real>& dst);
	static unsigned CalcLoopCount(const Config& cfg);
	double Freq2EnergyQuantile(double f);
	double EnergyQuantile2Freq(double q);
	void CreateDesign();
	void ReadDesign();
	void CreateTimeDomain(unique_fftw_arr<fftw_real>& sampbuf, double& norm, unique_fftwf_plan& plan);
	void CreatePCM(const unique_num_array<fftw_real>& sampbuf);
	virtual void Run();
	void DoSweep();
	/// Play OutBuf loopcount times or unless termrq.
	void PlayNoise(unsigned loopcount);
	void PlaySilence(unsigned loopcount);
	/// Move samples in OutBuf to the next channel.
	void NextChannel();
 public:
	/// Create output worker
	/// @param cfg global configuration.
	AnalyzeOut(const Config& cfg);
	/// Setup output worker.
	/// @return true if there is something to do asynchronously.
	bool Setup();
};

#endif // ANALYZE_H_
