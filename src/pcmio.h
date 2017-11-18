#ifndef PCMIO_H_
#define PCMIO_H_

#include "utils.h"

#include <stdlib.h>
#include <math.h>

enum class Format
{	I16
,	I16_SWAP
,	F32
};

typedef double fftw_real;

struct MinMax
{	double Min;
	double Max;
	void reset() { Min = +INFINITY; Max = -INFINITY; }
	double store(double val)
	{	if (val < Min)
			Min = val;
		if (val > Max)
			Max = val;
		return val;
	}
};

class PCMinput
{	PCMinput(const PCMinput&) = delete;
	const PCMinput& operator=(const PCMinput&) = delete;
 public: // public config
	const Format Fmt;
	const size_t BytesPerSample;
	const bool Differential;
	const size_t AddBins;
	const double Gain[2];
 private: // configuration
 	void (PCMinput::*const ReadFunc)();
 private: // working set
	size_t Len;
	const char* Sp;
	fftw_real* Dp[2];
	double Ch[2];
 public: // results
	MinMax Limits[2];
 private: // reader functions
	void LoadShort();
	void LoadShortSwap();
	void LoadFloat();
	void LoadShortAdd();
	void LoadShortSwapAdd();
	void LoadFloatAdd();
 private: // converter function
	void Store();
	void StoreAdd();
	void StoreDiff();
	void StoreDiffAdd();
 public: // public API
	PCMinput(Format fmt, bool differential = false, size_t binsize = 1, double (*gain)[2] = NULL);
	void convert(const scoped_array<fftw_real>& dst1, const scoped_array<fftw_real>& dst2, const scoped_array<char>& src, bool add);
	void reset();
	void discard(FILE* in, size_t count, const scoped_array<char>& buffer) const;
	void ASCIIdump(FILE* out, const scoped_array<char>& src) const;
};

class PCMoutput
{	const Format Fmt;
 public:
	PCMoutput(Format fmt) : Fmt(fmt) {}
};

#endif // PCMIO_H_
