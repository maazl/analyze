#ifndef PCMIO_H_
#define PCMIO_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cinttypes>

#include "unique_array.h"


struct Format
{	enum Enum
	{	I16
	,	I16_SWAP
	,	F32
	} Value;
	constexpr Format(Enum val) : Value(val) {}
	constexpr operator Enum() const { return Value; }
	constexpr size_t getSize() const { return Value == F32 ? sizeof(float) : sizeof(int16_t); }
	constexpr float getFSR() const { return Value == F32 ? 1. : 32768.; }
};

typedef float fftw_real;

struct MinMax
{	fftw_real Min;
	fftw_real Max;
	void reset() { Min = +INFINITY; Max = -INFINITY; }
	fftw_real store(fftw_real val)
	{	if (val < Min)
			Min = val;
		if (val > Max)
			Max = val;
		return val;
	}
};

class PCMIO
{	PCMIO(const PCMIO&) = delete;
	const PCMIO& operator=(const PCMIO&) = delete;
 public: // public config
	const Format Fmt;
	const size_t BytesPerSample;
 protected:
	PCMIO(Format fmt) : Fmt(fmt), BytesPerSample(2 * fmt.getSize()) {}
};

class PCMinput : public PCMIO
{public: // public config
	const bool Differential;
	const fftw_real Gain[2];
 private: // configuration
 	void (PCMinput::*const ReadFunc)();
 private: // working set
	size_t Len;
	const char* Sp;
	fftw_real* Dp[2];
	fftw_real Ch[2];
	unique_num_array<char> Buffer;
 public: // results
	MinMax Limits[2];
 private: // reader functions
	void LoadShort();
	void LoadShortSwap();
	void LoadFloat();
 private: // converter function
	void Store();
	void StoreAdd();
	void StoreDiff();
	void StoreDiffAdd();
 public: // public API
	PCMinput(Format fmt, bool differential = false, const fftw_real (*gain)[2] = NULL);
	void reset();
	void convert(const unique_num_array<fftw_real>& dst1, const unique_num_array<fftw_real>& dst2, const unique_num_array<char>& src, bool add = false);
	void discard(FILE* in, size_t count, const unique_num_array<char>& buffer) const;
	void discard(FILE* in, size_t count);
	void ASCIIdump(FILE* out, const unique_num_array<char>& src) const;
	void read(FILE* in, const unique_num_array<fftw_real>& dst1, const unique_num_array<fftw_real>& dst2, bool add = false, FILE* asciidump = NULL);
};

class PCMoutput : public PCMIO
{public: // public config
	const fftw_real Gain;
 private:
	unique_num_array<char> Buffer;
 public:
	PCMoutput(Format fmt, fftw_real gain = 1.);
	void convert(const unique_num_array<fftw_real>& src, const unique_num_array<char>& dst);
	void convert(const unique_num_array<fftw_real>& src1, const unique_num_array<fftw_real>& src2, const unique_num_array<char>& dst);
	void write(FILE* out, const unique_num_array<fftw_real>& src);
	void write(FILE* out, const unique_num_array<fftw_real>& srcL, const unique_num_array<fftw_real>& srcR);
	void zero(FILE* out, size_t nsamp);
	void WAVheader(FILE* out, size_t nsamp, size_t sfreq);
};

#endif // PCMIO_H_
