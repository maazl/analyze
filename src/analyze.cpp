#include "pca.h"
#include "parser.h"
#include "utils.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <stdarg.h>
#include <errno.h>

#include <fftw3.h>
#include <complex>
#include <vector>
using namespace std;
typedef double fftw_real;
typedef complex<double> Complex;

#define M_2PI (2.*M_PI)
#define M_3PI (3.*M_PI)
#ifndef M_PI_2
#define M_PI_2 (M_PI*.5)
#endif
#define M_PI_180 (M_PI/180.)
#define M_180_PI (180./M_PI)

#define N_MAX (65536*8)
#define CA_MAX 2
#define HA_MAX 5

// avoid name clash with math.h
#define fmin fmin__
#define fmax fmax__


// data buffers
static scoped_array<int16_t> inbuffertmp;      // Buffer for raw input
static scoped_fftw_arr<fftw_real> inbuffer[2]; // Buffer for nominator and denominator input
static scoped_fftw_arr<fftw_real> ovrbuffer[2];// Buffer for overridden nominator/denominator
static scoped_fftw_arr<fftw_real> outbuffer[2];// Buffer for FFT(inbuffer[])
static scoped_fftw_arr<fftw_real> ccbuffer1;   // Buffer for cross correlation temporary data
static scoped_fftw_arr<fftw_real> ccbuffer2;   // Buffer for cross correlation of outbuffer
static scoped_array<fftw_real> window;         // Buffer for window function
static scoped_array<int> harmonics;            // Buffer for harmonics dispatch table

static double noiselvl_;

static const double minval = 1E-20;

static double getweight(double a1, double a2, double)
{ /*a1 -= noiselvl;
 a2 -= noiselvl;
 if (a1 <= 0)
 a1 = 0;
 if (a2 <= 0)
 a2 = 0;*/
	return noiselvl_ / (1 / sqr(a1 + minval) + 1 / sqr(a2 + minval));
}

static double getweightD(double a1, double a2, double)
{ /*a1 -= noiselvl;
 a2 -= noiselvl;
 if (a1 <= 0)
 a1 = 0;
 if (a2 <= 0)
 a2 = 0;*/
	a1 += minval;
	double a1q = 1 / sqr(a1);
	return noiselvl_ / (a1q + sqrt(a1q + 1 / sqr(a1 + a2)) / (a2 + minval));
}

static double getconstweight(double, double, double)
{
	return noiselvl_;
}

static double get1_fweight(double, double, double f)
{
	return noiselvl_ / f;
}

static void vectorscale(const scoped_array<fftw_real>& data, double factor)
{	size_t len = data.size();
	fftw_real* dp = data.begin();
	while (len--)
		*dp++ *= factor;
}

static void vectoradd(const scoped_array<fftw_real>& dst, const scoped_array<fftw_real>& src)
{	size_t len = dst.size();
	assert(len == src.size());
	fftw_real* dp = dst.begin();
	const fftw_real* sp = src.begin();
	while (len--)
		*dp++ += *sp++;
}

static void vectormul(const scoped_array<fftw_real>& dst, const scoped_array<fftw_real>& src)
{	size_t len = dst.size();
	assert(len == src.size());
	fftw_real* dp = dst.begin();
	const fftw_real* sp = src.begin();
	while (len--)
		*dp++ *= *sp++;
}

static void vectordiv(const scoped_array<fftw_real>& dst, const scoped_array<fftw_real>& src)
{	size_t len = dst.size();
	assert(len == src.size());
	fftw_real* dp = dst.begin();
	const fftw_real* sp = src.begin();
	while (len--)
		*dp++ /= *sp++;
}

// config
static fftw_real gainadj[2] = { 1, 1 }; // gain {l, r}
static unsigned N = 8192;     // FFT length
static double noiselvl = 1;   // ?
static unsigned winfn = 0;    // window function: 0 = rectangle, 1 = Bartlett, 2 = Hanning, 3 = Hamming, 4 = Blackman, 5 = Blackman-Harris
static double freq = 48000;   // sampling rate
static double fmin = 0;       // minimum freuqncy for FFT analysis
static double fmax = INFINITY;// minimum freuqncy for FFT analysis
static double famin = 1;      // ignore frquencies below famin for calculation of screen output
static double famax = 1E99;   // ignore frquencies above famax for calculation of screen output
static bool writeraw = false; // write raw data to file
static bool writesrc = false; // write input data to file
static bool writedata = false;// write analysis data to file
static bool writewindow = false;// write window function data to file
static bool mpca = false;     // analysis method PCA
static bool mfft = false;     // analysis method FFT
static bool mxy = false;      // analysis method XY
static unsigned purgech = 1;  // set the first FFT frequencies to 0
static unsigned discsamp = 0; // skip the first samples
static bool disctrail = false;// comsume trailing samples after completion
static unsigned addch = 1;    // binsize in raw samples
static unsigned addloop = 1;  // add raw data before analysis # times
static bool incremental = false;// incremental mode (add all raw data)
static double rref = 1;       // value of the reference resistor in impedance measurements
static bool swapbytes = false;// swap bytes on PCM input
static unsigned scalemode = 1;// l/r matrix decoder: 1 = L=l & R=r, 2 = L=r & R=l-r, 3 = L=r & R=l
static bool stereo = false;   // Stereo aggregate mode (Toggle harmonics)
static unsigned loops = 1;    // number of analysis loops
static unsigned lpause = 10;  // number of loops between zero calibration parts
static unsigned zeromode = 0; // zero calibration mode: 0 = none, 1 = read, 2 = generate, 3 = generatedelta, 4 = generate part 2, 5 = generatedelta part 2
static unsigned gainmode = 0; // gain calibration mode: 0 = none, 1 = read, 2 = generate, 3 = generatedelta
static double linphase = 0;   // linear phase correction [s]
static bool normalize = false;// normalize L+R to 1. for impedance measurements
static unsigned binsz = 1;    // binsize in FFT channels
static double fbinsc = 0;     // logarithmic binsize: fmin/fmax = 1 + fbinsc
static double f_inc = 1;      // Absolute increment for harmonic table calculation
static double f_log = 0;      // Relative increment for harmonic table calculation
static unsigned harmonic = 0; // analyze up to # harmonics
static bool crosscorr = false;// Calculate and remove time delay by cross correlation
static const char* datafile = "data.dat"; // filename for analysis data
static const char* zerofile = "zero.dat"; // file name for zero calibration data
static const char* zerodifffile = "zeroD.dat"; // file name for differential zero calibration data
static const char* gainfile = "gain.dat"; // file name for gain calibration data
static const char* gaindifffile = "gainD.dat"; // file name for differential gain calibration data
static const char* rawfile = "raw.dat";  // file name for raw data
static const char* srcfile = "source.dat";// file name for input data
static const char* windowfile = "window.dat"; // file name for window data
static struct ovrwrt          // overwrite channel with ...
{	const char* file;           // ... file
	unsigned column;            // ... column in file
} overwrt[2] =
{	{ NULL, 1 }
,	{ NULL, 1 }
};
static const char* infile = "-";  // input file name
static const char* execcmd = NULL;// shell command to execute after analysis
static const char* plotcmd = NULL;// string to write to stdout after analysis
static double (*weightfn)(double, double, double) = getweight;// weight function
// internal vars
static fftw_plan P;    // FFT plan for forward transformation
static fftw_plan PI;   // FFT plan for inverse transformation

static void createwindow(const scoped_array<fftw_real>& dst, int type)
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
	vectorscale(dst, len / sum);
}

struct limits
{	int Min;
	int Max;
	void reset() { Min = INT_MAX; Max = INT_MIN; }
	int store(int val)
	{	if (val < Min)
			Min = val;
		if (val > Max)
			Max = val;
		return val;
	}
};
static limits minmax[2];

class reader16
{	size_t Len;
	const int16_t* Sp;
	fftw_real* Dp[2];
	double Ch[2];
	void read2float();
 public:
	reader16(const scoped_array<fftw_real>& dst1, const scoped_array<fftw_real>& dst2, const scoped_array<int16_t>& src);
	reader16(const reader16&) = delete;
	const reader16& operator=(const reader16&) = delete;
	void short2float2();
	void short2float2add();
	void short2float2window();
	void short2floatD();
	void short2floatDadd();
	void short2floatDwindow();
};

reader16::reader16(const scoped_array<fftw_real>& dst1, const scoped_array<fftw_real>& dst2, const scoped_array<int16_t>& src)
:	Len(dst1.size())
,	Sp(src.get())
,	Dp{dst1.get(), dst2.get()}
{	assert(2 * addch * Len == src.size() && Len == dst2.size());
}

void reader16::read2float()
{	Ch[0] = Ch[1] = 0;
	unsigned i = addch;
	if (swapbytes)
		do
		{	Ch[0] += minmax[0].store(bswap(Sp[0]));
			Ch[1] += minmax[1].store(bswap(Sp[1]));
			Sp += 2;
		} while (--i);
	else
		do
		{	Ch[0] += minmax[0].store(Sp[0]);
			Ch[1] += minmax[1].store(Sp[1]);
			Sp += 2;
		} while (--i);
}

void reader16::short2float2()
{	while (Len--)
	{	read2float();
		*Dp[0]++ = Ch[0] * gainadj[0];
		*Dp[1]++ = Ch[1] * gainadj[1];
	}
}

void reader16::short2float2add()
{	while (Len--)
	{	read2float();
		*Dp[0]++ += Ch[0] * gainadj[0];
		*Dp[1]++ += Ch[1] * gainadj[1];
	}
}

void reader16::short2float2window()
{	const fftw_real* win = window.get();
	while (Len--)
	{	read2float();
		*Dp[0]++ = Ch[0] * *win * gainadj[0];
		*Dp[1]++ = Ch[1] * *win++ * gainadj[1];
	}
}

void reader16::short2floatD()
{	while (Len--)
	{	read2float();
		*Dp[1]++ = Ch[0] * gainadj[0] - (*Dp[0]++ = Ch[1] * gainadj[1]);
	}
}

void reader16::short2floatDadd()
{	while (Len--)
	{	read2float();
		Ch[1] *= gainadj[1];
		*Dp[0]++ += Ch[1];
		*Dp[1]++ += Ch[0] * gainadj[0] - Ch[1];
	}
}

void reader16::short2floatDwindow()
{	const fftw_real* win = window.get();
	while (Len--)
	{	read2float();
		*Dp[1]++ = Ch[0] * *win * gainadj[0] - (*Dp[0]++ = Ch[1] * *win * gainadj[1]);
		++win;
	}
}


static void write1ch(FILE* out, const scoped_array<fftw_real>& data)
{	size_t len = data.size();
	const fftw_real* sp = data.get();
	while (len--)
		fprintf(out, "%g\n", *sp++);
}

static void write2ch(FILE* out, const scoped_array<int16_t>& data)
{	size_t len = data.size() >> 1;
	const int16_t* sp = data.get();
	if (swapbytes)
		while (len--)
		{	fprintf(out, "%i\t%i\n", bswap(sp[0]), bswap(sp[1]));
			sp += 2;
		}
	else
		while (len--)
		{	fprintf(out, "%i\t%i\n", sp[0], sp[1]);
			sp += 2;
		}
}
static void write2ch(FILE* out, const scoped_array<fftw_real>& data1, const scoped_array<fftw_real>& data2)
{	size_t len = data1.size();
	assert (len == data2.size());
	const fftw_real* sp1 = data1.get();
	const fftw_real* sp2 = data2.get();
	while (len--)
		fprintf(out, "%g\t%g\n", *sp1++, *sp2++);
}

/*static void writepolar(FILE* out, const fftw_real* data, size_t len, double inc)
{
	const fftw_real* data2 = data + len;
	// 1st line
	fprintf(out, "0\t%g\t0\n", *data++);
	len = 1;
	while (data < --data2)
		fprintf(out, "%g\t%g\t%g\n", len++ * inc, *data++, *data2);
}*/

static void writecomplex(FILE* out, const scoped_array<Complex>& data)
{	size_t len = data.size();
	const Complex* sp = data.get();
	while (len--)
	{	fprintf(out, "%14g\t%14g\t%14g\t%14g\n", sp->real(), sp->imag(), abs(*sp), arg(*sp) * M_180_PI);
		++sp;
	}
}

static void write4complex(FILE* out, const scoped_array<array<Complex,4>>& data)
{	size_t len = data.size();
	const array<Complex,4>* sp = data.get();
	while (len--)
	{	//Complex det = (*data)[0] * (*data)[3] - (*data)[1] * (*data)[2];
		fprintf(out, "%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\n",
			(*sp)[0].real(), (*sp)[0].imag(), (*sp)[1].real(), (*sp)[1].imag(),
			(*sp)[2].real(), (*sp)[2].imag(), (*sp)[3].real(), (*sp)[3].imag(),
			abs((*sp)[0]), arg((*sp)[0]) * M_180_PI, abs((*sp)[1]), arg((*sp)[1]) * M_180_PI,
			abs((*sp)[2]), arg((*sp)[2]) * M_180_PI, abs((*sp)[3]), arg((*sp)[3]) * M_180_PI);
		//, det.real(), det.imag());
		++sp;
	}
}

static void readcomplex(FILE* in, const scoped_array<Complex>& data)
{	size_t len = data.size();
	Complex* dp = data.get();
	while (len--)
	{	fscanf(in, "#%*[^\n]"); // skip comments
		double a, b;
		if (fscanf(in, "%lg%lg%*[^\n]", &a, &b) != 2)
			die(27, "Failed to read complex data: %s", strerror(errno));
		//(stderr, "%g\t%g\n", a,b);
		*dp++ = Complex(a, b);
	}
}

static void read4complex(FILE* in, const scoped_array<array<Complex,4>>& data)
{	size_t len = data.size();
	array<Complex,4>* dp = data.get();
	while (len--)
	{	fscanf(in, "#%*[^\n]"); // skip comments
		double a, b, c, d, e, f, g, h;
		if (fscanf(in, "%lg%lg%lg%lg%lg%lg%lg%lg%*[^\n]", &a, &b, &c, &d, &e, &f, &g, &h) != 8)
			die(27, "Failed to read four complex values: %s", strerror(errno));
		//(stderr, "%g\t%g\n", a,b);
		(*dp)[0] = Complex(a, b);
		(*dp)[1] = Complex(c, d);
		(*dp)[2] = Complex(e, f);
		(*dp)[3] = Complex(g, h);
		++dp;
	}
}

static void readfloat_2(FILE* in, unsigned column, const scoped_array<fftw_real>& dest)
{	size_t count = dest.size();
	fftw_real* dp = dest.get();
	while (count--)
	{	fscanf(in, "#%*[^\n]"); // skip comments
		unsigned col = column;
		while (--col)
			fscanf(in, "%*s");
		if (fscanf(in, "%lg%*[^\n]", dp) != 1)
			die(27, "Failed to read column %u from data file.", column);
		++dp;
	}
}

static void dofft()
{	// forwardtransformation
	fftw_execute_r2r(P, inbuffer[0].get(), outbuffer[0].get());
	fftw_execute_r2r(P, inbuffer[1].get(), outbuffer[1].get());

	static const double minscale = 1E-15;
	if (purgech)
	{	// purge DC (meaningless without DC-coupling)
		vectorscale(outbuffer[0].slice(0, purgech + 1), minscale);
		vectorscale(outbuffer[0].slice(N - purgech, purgech), minscale);
		vectorscale(outbuffer[1].slice(0, purgech + 1), minscale);
		vectorscale(outbuffer[1].slice(N - purgech, purgech), minscale);
	}
	// append constant zero to FFT result to simplify analysis
	outbuffer[0][N] = 0;
	outbuffer[1][N] = 0;

	// Calculate cross correlation to compensate for constant group delay.
	if (crosscorr)
	{  // calculate outbuffer1[] * conjugate(outbuffer2[])
	   // f(0)
		ccbuffer1[0] = outbuffer[0][0] * outbuffer[1][0];
		// f(1..N/2-1)
		for (unsigned i = 1; i < N / 2; ++i)
		{	Complex c = Complex(outbuffer[0][i], outbuffer[0][N - i]) * Complex(outbuffer[1][i], -outbuffer[1][N - i]);
			ccbuffer1[i] = c.real();
			ccbuffer1[N - i] = c.imag();
		}
		// f(N/2)
		ccbuffer1[N / 2] = outbuffer[0][N / 2] * outbuffer[1][N / 2];

		// do the cross correlation
		fftw_execute_r2r(PI, ccbuffer1.get(), ccbuffer2.get());
		/*FILE* F = fopen("cc.dat", "w");
		 write1ch(F, ccbuffer2, N);
		 fclose(F);*/

		// use the result
		double phiinc = M_2PI / N;
		double asum = 0;
		double bsum = 0;
		double sum = 0;
		for (unsigned i = 0; i < N; ++i)
		{
			double amp = ccbuffer2[i];
			amp *= amp;
			asum += cos(phiinc * i) * amp;
			bsum += sin(phiinc * i) * amp;
			sum += amp;
		}
		asum /= sum;
		bsum /= sum;
		/*fprintf(stderr, "***** %12g %12g %12g %12g\n",
		 sqrt(asum*asum+bsum*bsum), atan2(bsum, asum)*M_180_PI, asum, bsum);*/

		// calc linphase to compensate for the group delay
		linphase = atan2(bsum, asum) * N / freq;
	}
}

static scoped_array<double> wsums;
static scoped_array<Complex> gain;
static scoped_array<Complex> gainD;
static scoped_array<array<Complex,4>> zero;
static scoped_array<array<Complex,4>> zeroD;

// Calibration
static void docal(int bin, double f, Complex& U, Complex& I)
{
	switch (gainmode)
	{
	case 1: // read
		U /= gain[bin];
		/*t = I * gainD[len];
		 I -= U * gainD[len];
		 U -= t;*/
		break;
	case 2: // write
	case 3:
		{
			double weight = sqrt((*weightfn)(abs(U), abs(I), f));
			wsums[bin] += weight;
			U /= gain[bin];
			gainD[bin] += I / U * weight;
		}
	}
	if (normalize)
	{
		Complex s = 1. / (U + I);
		U *= s;
		I *= s;
	}
	if (zeromode & 1)
	{
		array<Complex,4>& cp = zero[bin];
		Complex t = U; // multiply (U,I) by (*cp)^(-1). The matrix inversion is easy because det(*cp) == 1.
		//Complex det = cp[0]*cp[3] - cp[1]*cp[2];
		U = (U * cp[3] - I * cp[1]);
		I = (-t * cp[2] + I * cp[0]);
		//U = (U * c22 - I * c12);
		//I = (-t * c21 + I * c11);
	}
	switch (zeromode)
	{
	case 2: // write part 1
	case 3:
		zeroD[bin][0] += U;
		zeroD[bin][2] += I;
		break;
	case 4: // write part 2
	case 5:
		zeroD[bin][1] += U;
		zeroD[bin][3] += I;
	}
}

// Phase unwrapper
// Adjusts phase by adding a multiple of 2pi so that the result
// is as close as possible to lph.
//   lph     last phase
//   phase   current calculated phase
//   returns unwrapped phase
static double unwrap(double lph, double phase)
{
	return !isfinite(lph) ? phase : phase - M_2PI * floor((phase - lph) / M_2PI + .5);
}

class FFTbin
{
public:
	enum StoreRet
	{	BelowMin,  // frequency less than fmin
		AboveMax,  // frequency above fmax
		Ready,     // calculated values available
		Aggregated,     // bin used for aggregation only
		Skip       // skip this bin because it is a harmonic
	};
private:
	struct aggentry
	{	double f;
		double Uabs; // Magnitude of nomiator
		double Uarg; // Phase of nominator
		double Iabs; // Magnitude of denominator
		double Iarg; // Phase of denomiator
		double Zabs; // Magnitude of quotient
		double Zarg; // Phase of quotient
		double D;    // Group delay
		double W;    // weight sum
		// internals
		double lf;   // last frequency (for numerical derivative)
		double lUarg;// last phase of nomiator
		double lIarg;// last phase of denomiator
		double lZarg;// last phase of quotient
		double fnext;// next frequency for bin size
		unsigned binc;// number of bins accumulated
	};
private:
	const double finc;

	aggentry agg[2 * HA_MAX + 1];
	int ch;
	aggentry* curagg;
	Complex Zcache;

public:
	FFTbin(double finc) : finc(finc)
	{	memset(agg, 0, sizeof agg);
	}
	StoreRet StoreBin(unsigned bin);

	static void PrintHdr(FILE* dst);
	void PrintBin(FILE* dst) const;

	double f() const    { return curagg->f; }    ///< frequency
	Complex U() const   { return polar(curagg->Uabs, curagg->Uarg); } ///< voltage, nominator or wanted signal
	double Uabs() const { return curagg->Uabs; } ///< voltage magnitude
	double Uarg() const { return curagg->Uarg; } ///< voltage phase
	Complex I() const   { return polar(curagg->Iabs, curagg->Iarg); } ///< current, denominator or reference signal
	double Iabs() const { return curagg->Iabs; } ///< current magnitude
	double Iarg() const { return curagg->Iarg; } ///< current phase
	Complex Z() const   { return Zcache; }       ///< impedance, quotient or relative signal
	double Zabs() const { return curagg->Zabs; } ///< impedance magnitude
	double Zarg() const { return curagg->Zarg; } ///< impedance phase
	double D() const    { return curagg->D / M_2PI; } ///< group delay
	double W() const    { return curagg->W; }    ///< weight
	int h() const       { return ch; }           ///< harmonic
};

FFTbin::StoreRet FFTbin::StoreBin(unsigned bin)
{
	curagg = NULL;
	// frequency
	double f = bin * finc;
	ch = harmonics[bin];
	if ((unsigned)(abs(ch) - 1) >= HA_MAX)
		return Skip;
	curagg = agg + ch + HA_MAX;
	unsigned base = bin;        //ch != 0 ? bin/abs(ch) : bin;
	// retrieve coefficients
	Complex U(outbuffer[0][bin], bin && bin != N / 2 ? outbuffer[0][N - bin] : 0);
	Complex I(outbuffer[1][base], bin && bin != N / 2 ? outbuffer[1][N - base] : 0);
	// calibration
	docal(bin, f, U, I);
	// phase correction
	U *= Complex(cos(linphase * f), sin(linphase * f));
	// calc Y
	Complex Z(U / I);
	// convert to polar
	double Uabs = abs(U);
	double Uarg = unwrap(curagg->lUarg, arg(U));
	double Iabs = abs(I);
	double Iarg = unwrap(curagg->lIarg, arg(I));
	double Zabs = abs(Z);
	double Zarg = unwrap(curagg->lZarg, arg(Z));
	// group delay
	double D = (Zarg - curagg->lZarg) / (f - curagg->lf);
	// store values for next point
	curagg->lUarg = Uarg;
	curagg->lIarg = Iarg;
	curagg->lZarg = Zarg;
	curagg->lf = f;

	// weight
	double w = (*weightfn)(Uabs, Iabs, f);
	if (curagg->binc == 0)
	{  // init
		curagg->f = f * w;
		curagg->Uabs = Uabs * w;
		curagg->Uarg = Uarg * w;
		curagg->Iabs = Iabs * w;
		curagg->Iarg = Iarg * w;
		curagg->Zabs = Zabs * w;
		curagg->Zarg = Zarg * w;
		curagg->D = D * w;
		curagg->W = w;
		curagg->fnext = f * (1 + fbinsc) - finc;
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

void FFTbin::PrintHdr(FILE* dst)
{
	fputs("#f\t|U|\targ U\t|I|\targ I\t|Z|\targ Z\tZ real\tZ imag\tweight\tdelay\tharmon.\n", dst);
}

void FFTbin::PrintBin(FILE* dst) const
{
	fprintf(dst, "%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%6i\n",
	// f    |Hl|      phil               |Hr|      phir
	   f(), Uabs(), Uarg() * M_180_PI, Iabs(), Iarg() * M_180_PI,
	// |Hl|/|Hr| phil-phir          re          im
	   Zabs(), Zarg() * M_180_PI, Z().real(), Z().imag(),
	// weight delay harmonic
	   W(), D(), h());
}

static const OptionDesc OptMap[] =
{	MkOpt("ainc", "incremental mode", &incremental, true)
,	MkOpt("al",   "average over multiple cycles of samples", &addloop)
,	MkOpt("bin",  "average FFT channels", &binsz)
,	MkOpt("bn",   "analysis block size", &N)
,	MkOpt("ca",   "add samples to bins", &addch)
,	MkOpt("df",   "name of the FFT data file", &datafile)
,	MkOpt("exec", "execute shell command after data available", &execcmd)
,	MkOpt("famax","upper frequency range for LCR analysis", &famax)
,	MkOpt("famin","lower frequency range for LCR analysis", &famin)
,	MkOpt("fbin", "average FFT channels with logarithmic bandwidth", &fbinsc)
,	MkOpt("finc", "linear increment for used FFT channels", &f_inc)
,	MkOpt("flog", "logarithmic increment for used FFT channels", &f_log)
,	MkOpt("fmax", "upper frequency range for analysis", &fmax)
,	MkOpt("fmin", "lower frequency range for analysis", &fmin)
,	MkOpt("fq",   "sampling frequency", &freq)
,	MkOpt("g2f",  "name of validation file of gain calibration", &gaindifffile)
,	MkOpt("gd",   "verify gain calibration", &gainmode, 3U)
,	MkOpt("gf",   "name of gain calibration file", &gainfile)
,	MkOpt("gg",   "generate gain calibration file", &gainmode, 2U)
,	MkOpt("gr",   "use gain calibration file", &gainmode, 1U)
,	MkOpt("h/f",  "use 1/f weight", &weightfn, get1_fweight)
,	MkOpt("harm", "take harmonics into account", &harmonic)
,	MkOpt("hd",   "use weight function for differential input mode", &weightfn, getweightD)
,	MkOpt("he",   "disable weight function", &weightfn, getconstweight)
,	MkOpt("in",   "name of input file", &infile)
,	MkOpt("ln",   "number of loops", &loops)
,	MkOpt("loop", "infinite number of loops", &loops, UINT_MAX)
,	MkOpt("lp",   "pause at matrix calibration", &lpause)
,	MkOpt("lvl",  "noise level for automatic weight function", &noiselvl)
,	MkOpt("mfft", "enable operation mode FFT", &mfft)
,	MkOpt("mpca", "enable operation mode PCA", &mpca)
,	MkOpt("mst",  "two channel mode", &stereo)
,	MkOpt("mxy",  "enable operation mode XY (preliminary)", &mxy)
,	MkOpt("olc",  "column to overwrite nominator", &overwrt[0].column)
,	MkOpt("olf",  "file name to overwrite nominator", &overwrt[0].file)
,	MkOpt("orc",  "column to overwrite denominator (reference)", &overwrt[1].column)
,	MkOpt("orf",  "file name to overwrite denominator (reference)", &overwrt[1].file)
,	MkOpt("pdc",  "purge first frequency channels", &purgech)
,	MkOpt("phcc", "fit group delay", &crosscorr)
,	MkOpt("phl",  "subtract constant group delay", &linphase)
,	MkOpt("plot", "write command to stdout after data available", &plotcmd)
,	MkOpt("psa",  "discard first samples", &discsamp)
,	MkOpt("pte",  "read input data till the end", &disctrail)
,	MkOpt("rf",   "name of raw data file", &rawfile)
,	MkOpt("rref", "reference resistor", &rref)
,	MkOpt("scm",  "input mode [0,2]", &scalemode, 0, 2)
,	MkOpt("sf",   "raw source data file name", &srcfile)
,	MkOpt("wd",   "(over)write FFT data file on the fly", &writedata)
,	MkOpt("wf",   "name of window function file", &windowfile)
,	MkOpt("win",  "select window function [0,5]", &winfn, 0, 5)
,	MkOpt("wr",   "write raw data", &writeraw)
,	MkOpt("ws",   "write source data file", &writesrc)
,	MkOpt("ww",   "write window function", &writewindow)
,	MkOpt("xb" ,  "swap bytes", &swapbytes)
,	MkOpt("z2f",  "name of validation file of matrix calibration", &zerodifffile)
,	MkOpt("zd",   "validate matrix calibration", &zeromode, 3U)
,	MkOpt("zf",   "name of matrix calibration file", &zerofile)
,	MkOpt("zg",   "generate matrix calibration file", &zeromode, 2U)
,	MkOpt("zn",   "normalize amplitudes", &normalize)
,	MkOpt("zr",   "use matrix calibration file", &zeromode, 1U)
};

int main(int argc, char* argv[])
{
	/*for (int i = 0; i < 3; ++i)
	 { puts("");
	 fflush(stdout);
	 DosSleep(500);
	 }*/

	// parse cmdl
	{	Parser parser(OptMap);
		while (--argc)
			parser.HandleArg(*++argv);
	}

	if (N > N_MAX)
		die(32, "FFT Length too large.");
	if (N * addch > N_MAX * CA_MAX)
		die(32, "Input data length too large.");
	if (overwrt[0].file && overwrt[1].file)
		infile = NULL;
	if (mxy & (mpca|mfft))
		die(34, "Invalid combination of measurement modes, e.g. FFT and XY.");

	// prepare some global vars
	noiselvl_ = 1 / noiselvl;
	freq /= addch;
	linphase *= M_2PI;
	f_inc -= .5;
	f_log += 1;
	gainadj[0] /= 32767;
	gainadj[1] /= 32767;
	// allocate buffers
	inbuffertmp.reset(2 * N * addch);
	window.reset(N);
	if (mxy)
	{	// reserve space for integrals and differentials too
		inbuffer[0].reset(3 * N);
		inbuffer[1].reset(3 * N);
		outbuffer[0].reset(3 * N + 1);
		outbuffer[1].reset(3 * N + 1);
	} else
	{	inbuffer[0].reset(N);
		inbuffer[1].reset(N);
		outbuffer[0].reset(N + 1);
		outbuffer[1].reset(N + 1);
	}
	if (crosscorr && mpca)
	{	ccbuffer1.reset(N);
		ccbuffer2.reset(N);
	}
	harmonics.reset(N / 2 + 1);
	wsums.reset(N / 2 + 1);
	gain.reset(N / 2 + 1);
	gainD.reset(N / 2 + 1);
	zero.reset(N / 2 + 1);
	zeroD.reset(N / 2 + 1);
	// create plan
	// fftw_real in[N], tout[N], power_spectrum[N/2+1];
	P = fftw_plan_r2r_1d(N, inbuffer[0].get(), outbuffer[0].get(), FFTW_R2HC, FFTW_ESTIMATE);
	PI = fftw_plan_r2r_1d(N, inbuffer[0].get(), outbuffer[0].get(), FFTW_HC2R, FFTW_ESTIMATE);
	// adjust fmax
	if (fmax > freq/2)
		fmax = freq/2;
	// create harmonics table
	{	harmonics.clear();
		int sign = 1;
		for (unsigned i = (int)floor(fmin / freq * N + .5); i <= floor(fmax / freq * N + .5); ++i)
		{	if (i)
			{
				for (unsigned j = 1; j <= harmonic && i * j <= N / 2; ++j)
					if (harmonics[i * j])
						goto next_f;
				for (unsigned j = 1; i * j <= N / 2; ++j)
					harmonics[i * j] = j * sign;
			}
			if (stereo)
				sign = -sign;
			i = (int)floor(i * f_log + f_inc);
		 next_f: ;
		}
		/*FILE* f = fopen("harm.dat", "w");
		 for (unsigned i = 0; i <= N/2; ++i)
		 fprintf(f, "%12g %8i\n", i*freq/N, harmonics[i]);
		 fclose(f);*/
	}

	createwindow(window, winfn);
	// write window data
	if (writewindow)
		write1ch(FILEguard(windowfile, "wt"), window);

	FILEguard in = NULL;
	if (infile)
	{	in = strcmp(infile, "-") == 0
			? binmode(stdin) // streaming
			: checkedopen(infile, "rb");
		// discard first samples
		while (discsamp > inbuffertmp.size())
		{	fread2(inbuffertmp.get(), 2 * sizeof inbuffertmp[0], inbuffertmp.size(), in);
			discsamp -= inbuffertmp.size();
		}
		fread2(inbuffertmp.get(), 2 * sizeof inbuffertmp[0], discsamp, in);
	}

	if (overwrt[0].file)
	{	ovrbuffer[0].reset(N);
		readfloat_2(FILEguard(overwrt[0].file, "r"), overwrt[0].column, ovrbuffer[0]);
	}
	if (overwrt[1].file)
	{	ovrbuffer[1].reset(N);
		readfloat_2(FILEguard(overwrt[1].file, "r"), overwrt[0].column, ovrbuffer[1]);
	}

	const scoped_array<fftw_real> input[2] =
	{	scoped_array<fftw_real>(inbuffer[0].slice(0, N))
	,	scoped_array<fftw_real>(inbuffer[1].slice(0, N))
	};

	input[0].clear(); // init with 0 because of incremental mode
	input[1].clear();
	wsums.clear();
	gainD.clear();
	zeroD.clear();
	// prepare gainmode
	switch (gainmode)
	{case 1: // read
	 case 3:
		readcomplex(FILEguard(gainfile, "r"), gain);
	}

 restart_zero:
	// prepare zeromode
	switch (zeromode)
	{case 1: // read
	 case 3:
		read4complex(FILEguard(zerofile, "r"), zero);
	}

	// operation loop
	unsigned loop = loops;
	unsigned addloops = 0;
	do
	{
		if (in)
		{
			fread2(inbuffertmp.get(), 2 * sizeof inbuffertmp[0] * addch, N, in);

			// write raw data
			if (writeraw)
				write2ch(FILEguard(rawfile, "wt"), inbuffertmp);

			// reset min/max
			minmax[0].reset();
			minmax[1].reset();

			reader16 rdr(input[scalemode >= 3], input[scalemode < 3], inbuffertmp);
			if (scalemode & 1)
			{	if (incremental || addloops)
					rdr.short2float2add();
				else if (winfn && addloop == 1)
					rdr.short2float2window();
				else
					rdr.short2float2();
			} else
			{	if (incremental || addloops)
					rdr.short2floatDadd();
				else if (winfn && addloop == 1)
					rdr.short2floatDwindow();
				else
					rdr.short2floatD();
			}
			// write raw status
			fprintf(stderr, "\nmin:\t%i\t%i\nmax:\t%i\t%i\n", minmax[0].Min, minmax[1].Min, minmax[0].Max, minmax[1].Max);
		}

		if (incremental || addloops)
		{	if (ovrbuffer[0])
				vectoradd(input[0], ovrbuffer[0]);
			if (ovrbuffer[1])
				vectoradd(input[1], ovrbuffer[1]);
		} else
		{	if (ovrbuffer[0])
				input[0].copyfrom(ovrbuffer[0]);
			if (ovrbuffer[1])
				input[1].copyfrom(ovrbuffer[1]);
		}

		// write source data
		if (writesrc)
			write2ch(FILEguard(srcfile, "wt"), inbuffer[0], inbuffer[1]);

		if (++addloops < addloop)
			continue; // add more data
		if (winfn && (incremental || addloops))
		{	// optimization: apply window function later
			vectormul(input[0], window);
			vectormul(input[1], window);
		}
		addloops = 0;

		// write raw data
		/*if (writeraw)
		 {  tout = fopen(rawfile, "wt");
		 write2ch(tout, inbuffer1, inbuffer2, N);
		 fclose(tout);
		 }*/

		if (mfft & mpca)
		{	// FFT
			dofft();

			vectorscale(outbuffer[0], sqrt(1. / N) / addch);
			vectorscale(outbuffer[1], sqrt(1. / N) / addch);

			// write data
			FILEguard tout = NULL;
			if (writedata)
			{	tout = checkedopen(datafile, "wt");
				FFTbin::PrintHdr(tout);
			}

			PCA<2> pcaRe;
			PCA<3> pcaIm;
			double PCAdataRe[2];
			double PCAdataIm[3];
			// some values are const
			PCAdataRe[1] = 1;
			//PCAdataIm[1] = 1;

			FFTbin calc(freq / N);

			// 1st line
			for (size_t len = 0; len <= N / 2; ++len)
			{  // do calculations and aggregations
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

				if (calc.f() / calc.h() < famin && calc.f() / calc.h() >= famax)
					continue;
				// component analysis
				PCAdataRe[0] = calc.Z().real();
				//PCAdataRe[2] = 1/af;
				//PCAdataRe[3] = f;
				PCAdataIm[0] = calc.Z().imag(); // fit imaginary part in conductivity
				PCAdataIm[1] = 1 / calc.f();
				PCAdataIm[2] = calc.f();
				//PCAdataIm[3] = 1/af;
				//printf("Re: %12g %12g %12g %12g\n", PCAdataRe[0], PCAdataRe[1], PCAdataRe[2], weight);

				// add values
				pcaRe.Store(PCAdataRe, calc.W());
				pcaIm.Store(PCAdataIm, calc.W());
			}

			// calculate summary
			Vektor<1> resRe = pcaRe.Result();
			Vektor<2> resIm = pcaIm.Result();

			fprintf(stderr, "resRe: %12g %12g %12g\n", resRe[0], resRe[1], resRe[2]);
			//printf("resRe: %12g %12g %12g %12g\n", resRe[0], resRe[1], resRe[2], resRe[3]);
			fprintf(stderr, "resIm: %12g %12g %12g %12g\n", resIm[0], resIm[1], resIm[2], resIm[3]);

			double R0 = rref * resRe[0];
			//double R1_f = rref*resRe[1];
			double C0 = -1 / M_2PI / resIm[0] / rref;
			double L0 = resIm[1] / M_2PI * rref;
			//double C1f = -resIm[0] / M_2PI / rref;

			// write summary
			fprintf(stderr, "\nreal: R [Ohm]     \t%12g\n"
			//       "      R/f [Ohm/Hz]\t%12g\t%12g @100Hz\n"
			    "imaginary: C [�F] \t%12g\n"
			//  "      Cf [F Hz]   \t%12g\t%12g @100Hz\n"
			    "imaginary: L [�H] \t%12g\n"
			//             , R0, R1_f, R1_f/100, C0, C1f, C1f*100);
			, R0, C0 * 1E6, L0 * 1E6);
			if (crosscorr)
				fprintf(stderr, "delay        \t%12g\n", linphase / M_2PI);
		}
		else if (mpca)
		{
			PCA<5> pca;
			double data[6];
			fftw_real* U = inbuffer[0].get() + 1;
			fftw_real* I = inbuffer[2].get() + 1;
			data[2] = 1;   // offset
			data[3] = 0;   // integral
			data[4] = 0;   // linear
			data[5] = 0;   // differential
			unsigned i = N - 3;
			do
			{
				data[0] = U[0] + U[1];
				data[1] = I[0] + I[1];
				data[3] += I[-1] + I[0];
				data[5] = I[-1] + I[0] - I[1] - I[2];
				pca.Store(*(double (*)[5])&data, (*weightfn)(data[0], data[1], i));
				data[4]++;
				U += 2;
				I += 2;
				i -= 2;
			} while (i > 0);

			Vektor<4> res = pca.Result() * rref;
			fprintf(stderr, "\nPCA: %12g %12g %12g %12g %12g %12g\n",
				res[0], res[1], 2. / freq / res[2], res[3], 1. / 2 * freq * res[2] / M_2PI / res[0], freq * res[4] / 2);
		}
		else if (mfft)
		{
			dofft();

			vectorscale(outbuffer[0].slice(0, N), sqrt(1. / N) / addch);
			vectorscale(outbuffer[1].slice(0, N), sqrt(1. / N) / addch);

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
			if (writedata)
			{	tout = checkedopen(datafile, "wt");
				FFTbin::PrintHdr(tout);
			}

			FFTbin calc(freq / N);

			for (size_t len = 0; len <= N / 2; ++len)
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

				if (calc.f() / calc.h() < famin && calc.f() / calc.h() >= famax)
					continue;
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
			fprintf(stderr, "\nreal (R)     \t%12g � %g\n"
					"imaginary (C)\t%12g � %g\n"
					"imaginary (L)\t%12g � %g\n", rref * R, rref * RE, 1 / (rref * C * M_2PI), 1 / (CE * rref * M_2PI), rref * L / M_2PI, LE * rref / M_2PI);
			if (crosscorr)
				fprintf(stderr, "delay        \t%12g\n", linphase / M_2PI);
		}
		else if (mxy)
		{	// we need to do an FFT here, at least for the calibration
			dofft();
			// normalize (including retransformation)
			vectorscale(outbuffer[0].slice(0, N), sqrt(1. / N / N) / addch);
			vectorscale(outbuffer[1].slice(0, N), sqrt(1. / N / N) / addch);

			const double inc = freq / N;
			// U(f)
			fftw_real* a1 = outbuffer[0].get();
			fftw_real* a2 = outbuffer[1].get();
			fftw_real* b1 = a1 + N;
			fftw_real* b2 = a2 + N;
			*b1 = 0; // well, somewhat easier this way
			*b2 = 0;
			for (int len = 0; a1 < b1; len += harmonic, a1 += harmonic, a2 += harmonic, b1 -= harmonic, b2 -= harmonic)
			{  // calc
				double f = len * inc;
				Complex U(*a1, *b1);
				Complex I(*a2, *b2);
				// calibration
				docal(len, f, U, I);
				// store data
				*a1 = U.real();
				*b1 = U.imag();
				*a2 = I.real();
				*b2 = I.imag();
				// The integrals and differentials are calculated in the frequency domain.
				// This is at the expence of approximate a factor 2 computing time, since we have to do
				// one forward and one backward transformation for the zero compensation anyway.
				// The advantage is that we do not have to deal with the 1/2 time slot linear phse shift
				// of the numeric integration/differentiation in the time domain.
				Complex di(0, f); // differential operator
				// store integrals
				Complex C = U / di;
				a1[N] = C.real();
				b1[N] = C.imag();
				C = I / di;
				a2[N] = C.real();
				b2[N] = C.imag();
				// store differentials
				C = U * di;
				a1[2 * N] = C.real();
				b1[2 * N] = C.imag();
				C = I * di;
				a2[2 * N] = C.real();
				b2[2 * N] = C.imag();
			}
			// clear DC and Nyquist frequencies of integral and differential, since they do not allow 90° phase shift.
			outbuffer[0][N] = outbuffer[1][N] = 0;
			outbuffer[0][2 * N] = outbuffer[1][2 * N] = 0;
			outbuffer[0][N + N / 2] = outbuffer[1][N + N / 2] = 0;
			outbuffer[0][2 * N + N / 2] = outbuffer[1][2 * N + N / 2] = 0;

			// now do the inverse transform to get the corrected data back.
			fftw_execute_r2r(PI, outbuffer[0].get(), inbuffer[0].get());
			fftw_execute_r2r(PI, outbuffer[0].get() + N, inbuffer[0].get() + N);
			fftw_execute_r2r(PI, outbuffer[0].get() + 2 * N, inbuffer[0].get() + 2 * N);
			fftw_execute_r2r(PI, outbuffer[1].get(), inbuffer[1].get());
			fftw_execute_r2r(PI, outbuffer[1].get() + N, inbuffer[1].get() + N);
			fftw_execute_r2r(PI, outbuffer[1].get() + 2 * N, inbuffer[1].get() + 2 * N);

			// write data
			if (writedata)
			{	FILEguard fout(datafile, "wt");
				const fftw_real* Up = inbuffer[0].get();
				const fftw_real* Ip = inbuffer[1].get();
				fputs("#t\tU\tI\t∫ U\t∫ I\tΔ U\t ΔI\n", fout);
				for (unsigned len = 0; len < N; ++len, ++Up, ++Ip)
					fprintf(fout, "%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\n",
						len / freq * harmonic, *Up, *Ip, Up[N], Ip[N], Up[2 * N], Ip[2 * N]);
			}
		}

		if (execcmd)
			system(execcmd);
		if (plotcmd)
		{	// for gnuplot!
			puts(plotcmd);
			fflush(stdout);
			/*#ifdef __OS2__
			 // Fix for curious bug in OS/2 preventing the fflush in the pipe from working reliable.
			 puts("\r\n");
			 fflush(stdout);
			 DosSleep(500);
			 #endif*/
		}

		// undo window function because of incremental mode
		// TODO: this causes a loss of precision!
		if (winfn && incremental)
		{	vectordiv(input[0], window);
			vectordiv(input[1], window);
		}

	} while (--loop);

	if (disctrail)
		fread2(inbuffertmp.get(), 2 * sizeof inbuffertmp[0] * addch, N, in);

	switch (gainmode)
	{case 2: // write
	 case 3:
		double* wp = wsums.end();
		for (Complex* cp = gainD.end(); --cp >= gainD.begin();)
			*cp /= *--wp;  // scale 2 average
		FILEguard fz(gainmode == 3 ? gaindifffile : gainfile, "wb");
		fputs("#f\treal\timag\tabs\targ\n", fz);
		writecomplex(fz, gainD);
	}
	switch (zeromode)
	{case 2: // write
	 case 3:
		zeromode += 2;
		puts("Zeromode calibration part one is now complete.\n"
				"Part 2 will start at the end of the contdown.\7");
		for (int loop = lpause; loop;)
		{	fprintf(stderr, "\r%u ", loop);
			fread2(inbuffertmp.get(), 2 * sizeof inbuffertmp[0] * addch, N, in);
			--loop;
		}
		puts("\nNow at part 2.");
		goto restart_zero;
	 case 4: // part 2
	 case 5:
		for (array<Complex,4>* cp = zeroD.end(); --cp >= zeroD.begin();)
		{	// scale to fit det *cp == 1
			Complex det = 1. / sqrt((*cp)[0] * (*cp)[3] - (*cp)[1] * (*cp)[2]);
			if ((((*cp)[0] + (*cp)[3]) * det).real() < 0)
				det = -det;
			(*cp)[0] *= det;
			(*cp)[1] *= det;
			(*cp)[2] *= det;
			(*cp)[3] *= det;
		}
		FILEguard fz(zeromode == 5 ? zerodifffile : zerofile, "wb");
		fputs("#f\tU->U re\tU->U im\tU->I re\tU->I im\tI->U re\tI->U im\tI->I re\tI->I im\n", fz);
		write4complex(fz, zeroD);
	}

	puts("completed.");
	// read until EOF
	if (disctrail)
		while (fread(inbuffertmp.get(), sizeof inbuffertmp[0], inbuffertmp.size(), in) > 0);

	return 0;
}

