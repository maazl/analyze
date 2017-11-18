#include <inttypes.h>
#include "pcmio.h"

using namespace std;

PCMinput::PCMinput(Format fmt, bool differential, size_t binsize, double (*gain)[2])
:	Fmt(fmt)
,	BytesPerSample(fmt == Format::F32 ? 2*sizeof(float) : 2*sizeof(int16_t))
,	Differential(differential)
,	AddBins(binsize)
,	Gain{gain ? (*gain)[0] : 1, gain ? (*gain)[1] : 1}
,	ReadFunc(array<void (PCMinput::*)(),6>{ &PCMinput::LoadShort, &PCMinput::LoadShortSwap, &PCMinput::LoadFloat, &PCMinput::LoadShortAdd, &PCMinput::LoadShortSwapAdd, &PCMinput::LoadFloatAdd }[(int)fmt + 3*(binsize>1)])
{	reset();
}

void PCMinput::convert(const scoped_array<fftw_real>& dst1, const scoped_array<fftw_real>& dst2, const scoped_array<char>& src, bool add)
{	Len = dst1.size();
	assert(BytesPerSample * AddBins * Len == src.size() && Len == dst2.size());
	Sp = src.get();
	Dp[0] = dst1.get();
	Dp[1] = dst2.get();
	if (Differential)
	{	if (add)
			StoreDiffAdd();
		else
			StoreDiff();
	} else
	{	if (add)
			StoreAdd();
		else
			Store();
	}
}

void PCMinput::reset()
{	Limits[0].reset();
	Limits[1].reset();
}

void PCMinput::LoadShort()
{	Ch[0] = Limits[0].store(((const int16_t*)Sp)[0]);
	Ch[1] = Limits[1].store(((const int16_t*)Sp)[1]);
	Sp += 2 * sizeof(int16_t);
}
void PCMinput::LoadShortAdd()
{	Ch[0] = Ch[1] = 0;
	unsigned i = AddBins;
	do
	{	Ch[0] += Limits[0].store(((const int16_t*)Sp)[0]);
		Ch[1] += Limits[1].store(((const int16_t*)Sp)[1]);
		Sp += 2 * sizeof(int16_t);
	} while (--i);
}
void PCMinput::LoadShortSwap()
{	Ch[0] = Limits[0].store(bswap(((const int16_t*)Sp)[0]));
	Ch[1] = Limits[1].store(bswap(((const int16_t*)Sp)[1]));
	Sp += 2 * sizeof(int16_t);
}
void PCMinput::LoadShortSwapAdd()
{	Ch[0] = Ch[1] = 0;
	unsigned i = AddBins;
	do
	{	Ch[0] += Limits[0].store(bswap(((const int16_t*)Sp)[0]));
		Ch[1] += Limits[1].store(bswap(((const int16_t*)Sp)[1]));
		Sp += 2 * sizeof(int16_t);
	} while (--i);
}
void PCMinput::LoadFloat()
{	Ch[0] = Limits[0].store(((const float*)Sp)[0]);
	Ch[1] = Limits[1].store(((const float*)Sp)[1]);
	Sp += 2 * sizeof(float);
}
void PCMinput::LoadFloatAdd()
{	Ch[0] = Ch[1] = 0;
	unsigned i = AddBins;
	do
	{	Ch[0] += Limits[0].store(bswap(((const float*)Sp)[0]));
		Ch[1] += Limits[1].store(bswap(((const float*)Sp)[1]));
		Sp += 2 * sizeof(float);
	} while (--i);
}

void PCMinput::Store()
{	while (Len--)
	{	(this->*ReadFunc)();
		*Dp[0]++ = Ch[0] * Gain[0];
		*Dp[1]++ = Ch[1] * Gain[1];
	}
}
void PCMinput::StoreAdd()
{	while (Len--)
	{	(this->*ReadFunc)();
		*Dp[0]++ += Ch[0] * Gain[0];
		*Dp[1]++ += Ch[1] * Gain[1];
	}
}
void PCMinput::StoreDiff()
{	while (Len--)
	{	(this->*ReadFunc)();
		*Dp[1]++ = Ch[1] * Gain[1] - (*Dp[0]++ = Ch[0] * Gain[0]);
	}
}
void PCMinput::StoreDiffAdd()
{	while (Len--)
	{	(this->*ReadFunc)();
		double ch0 = Ch[0] * Gain[0];
		*Dp[0]++ += ch0;
		*Dp[1]++ += Ch[1] * Gain[1] - ch0;
	}
}

void PCMinput::discard(FILE* in, size_t count, const scoped_array<char>& buffer) const
{	count *= BytesPerSample;
	while (count > buffer.size())
	{	fread2(buffer.get(), buffer.size(), 1, in);
		count -= buffer.size();
	}
	fread2(buffer.get(), count, 1, in);
}

void PCMinput::ASCIIdump(FILE* out, const scoped_array<char>& src) const
{	if (Fmt == Format::F32)
	{	size_t len = src.size() / sizeof(float);
		const float* sp = (float*)src.get();
		while (len--)
		{	fprintf(out, "%g\t%g\n", sp[0], sp[1]);
			sp += 2;
		}
	} else
	{	size_t len = src.size() / sizeof(int16_t);
		const int16_t* sp = (int16_t*)src.get();
		if (Fmt == Format::I16_SWAP)
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
}
