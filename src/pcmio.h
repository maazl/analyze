#ifndef PCMIO_H_
#define PCMIO_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cinttypes>

#include "unique_array.h"


enum struct Format : char
{	I16 = 4
,	I16_SWAP
,	I24
,	I24_SWAP
,	I32
,	I32_SWAP
,	F32
,	F32_SWAP
};
constexpr size_t getSize(Format value) { return value >= Format::F32 ? 4 : (size_t)value >> 1; }
constexpr float getFSR(Format value) { return value >= Format::F32 ? 1.F : ldexp(.5F, ((int)value & -2) << 2); }

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
	union endian_conv
	{	int32_t i;
		char    c[4];
	};
	PCMIO(Format fmt) : Fmt(fmt), BytesPerSample(2 * getSize(fmt)) {}
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
	void Load24();
	void Load24Swap();
	void Load24B();
	void Load24SwapB();
	void Load32();
	void Load32Swap();
	void LoadFloat();
	void LoadFloatSwap();
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
	void ASCIIdump(FILE* out, const unique_num_array<char>& src);
	void read(FILE* in, const unique_num_array<fftw_real>& dst1, const unique_num_array<fftw_real>& dst2, bool add = false, FILE* asciidump = NULL);
};

class PCMoutput : public PCMIO
{public: // public config
	const fftw_real Gain;
	const signed char SignCh2;
 private:
	unique_num_array<char> Buffer;
 public:
	/// Setup PCM converter
	/// @param fmt Sample format
	/// @param gain Gain factor
	/// @param signch2 Sign of channel 2 samples: 1 = same as channel 1, -1 = inverse, 0 = channel 1 only
	PCMoutput(Format fmt, fftw_real gain = 1., signed char signch2 = 1);
	void convert(const unique_num_array<fftw_real>& src, const unique_num_array<char>& dst);
	void convert(const unique_num_array<fftw_real>& src1, const unique_num_array<fftw_real>& src2, const unique_num_array<char>& dst);
	void write(FILE* out, const unique_num_array<fftw_real>& src);
	void write(FILE* out, const unique_num_array<fftw_real>& srcL, const unique_num_array<fftw_real>& srcR);
	void zero(FILE* out, size_t nsamp);
	void WAVheader(FILE* out, size_t nsamp, size_t sfreq);
};

#endif // PCMIO_H_
