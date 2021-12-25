#include "pcmio.h"
#include "utils.h"
#include <assert.h>
#include <stdio.h>

using namespace std;

#define BUFFERSIZE 65536


PCMinput::PCMinput(Format fmt, bool differential, const fftw_real (*gain)[2])
:	PCMIO(fmt)
,	Differential(differential)
,	Gain{ (gain ? (*gain)[0] : 1) / getFSR(fmt), (gain ? (*gain)[1] : 1) / getFSR(fmt) }
,	ReadFunc((is_big_endian()
	?	array<void (PCMinput::*)(),8>
		{	&PCMinput::LoadShort
		,	&PCMinput::LoadShortSwap
		,	&PCMinput::Load24B
		,	&PCMinput::Load24SwapB
		,	&PCMinput::Load32
		,	&PCMinput::Load32Swap
		,	&PCMinput::LoadFloat
		,	&PCMinput::LoadFloatSwap
		}
	:	array<void (PCMinput::*)(),8>
		{	&PCMinput::LoadShort
		,	&PCMinput::LoadShortSwap
		,	&PCMinput::Load24
		,	&PCMinput::Load24Swap
		,	&PCMinput::Load32
		,	&PCMinput::Load32Swap
		,	&PCMinput::LoadFloat
		,	&PCMinput::LoadFloatSwap
		})[(int)fmt - 4])
{	reset();
}

void PCMinput::reset()
{	Limits[0].reset();
	Limits[1].reset();
}

void PCMinput::convert(const unique_num_array<fftw_real>& dst1, const unique_num_array<fftw_real>& dst2, const unique_num_array<char>& src, bool add)
{	Len = dst1.size();
	assert(BytesPerSample * Len == src.size() && Len == dst2.size());
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

void PCMinput::LoadShort()
{	Ch[0] = Limits[0].store(((const int16_t*)Sp)[0]);
	Ch[1] = Limits[1].store(((const int16_t*)Sp)[1]);
	Sp += 2 * sizeof(int16_t);
}
void PCMinput::LoadShortSwap()
{	Ch[0] = Limits[0].store(bswap(((const int16_t*)Sp)[0]));
	Ch[1] = Limits[1].store(bswap(((const int16_t*)Sp)[1]));
	Sp += 2 * sizeof(int16_t);
}
void PCMinput::Load24()
{	endian_conv e;
	e.c[0] = Sp[0];
	e.c[1] = Sp[1];
	e.c[2] = Sp[2];
	e.c[3] = 0;
	Ch[0] = Limits[0].store(e.i);
	e.c[0] = Sp[3];
	e.c[1] = Sp[4];
	e.c[2] = Sp[5];
	Ch[1] = Limits[1].store(e.i);
	Sp += 6;
}
void PCMinput::Load24Swap()
{	endian_conv e;
	e.c[2] = Sp[0];
	e.c[1] = Sp[1];
	e.c[0] = Sp[2];
	e.c[3] = 0;
	Ch[0] = Limits[0].store(e.i);
	e.c[2] = Sp[3];
	e.c[1] = Sp[4];
	e.c[0] = Sp[5];
	Ch[1] = Limits[1].store(e.i);
	Sp += 6;
}
void PCMinput::Load24B()
{	endian_conv e;
	e.c[0] = 0;
	e.c[1] = Sp[0];
	e.c[2] = Sp[1];
	e.c[3] = Sp[2];
	Ch[0] = Limits[0].store(e.i);
	e.c[1] = Sp[3];
	e.c[2] = Sp[4];
	e.c[3] = Sp[5];
	Ch[1] = Limits[1].store(e.i);
	Sp += 6;
}
void PCMinput::Load24SwapB()
{	endian_conv e;
	e.c[0] = 0;
	e.c[3] = Sp[0];
	e.c[2] = Sp[1];
	e.c[1] = Sp[2];
	Ch[0] = Limits[0].store(e.i);
	e.c[3] = Sp[3];
	e.c[2] = Sp[4];
	e.c[1] = Sp[5];
	Ch[1] = Limits[1].store(e.i);
	Sp += 6;
}
void PCMinput::Load32()
{	Ch[0] = Limits[0].store(((int32_t*)Sp)[0]);
	Ch[1] = Limits[1].store(((int32_t*)Sp)[1]);
	Sp += 2 * sizeof(int32_t);
}
void PCMinput::Load32Swap()
{	int32_t f = ((const int32_t*)Sp)[0];
	bswap(&f);
	Ch[0] = Limits[0].store(f);
	f = ((const float*)Sp)[1];
	bswap(&f);
	Ch[1] = Limits[1].store(f);
	Sp += 2 * sizeof(float);
}
void PCMinput::LoadFloat()
{	Ch[0] = Limits[0].store(((const float*)Sp)[0]);
	Ch[1] = Limits[1].store(((const float*)Sp)[1]);
	Sp += 2 * sizeof(float);
}
void PCMinput::LoadFloatSwap()
{	float f = ((const float*)Sp)[0];
	bswap(&f);
	Ch[0] = Limits[0].store(f);
	f = ((const float*)Sp)[1];
	bswap(&f);
	Ch[1] = Limits[1].store(f);
	Sp += 2 * sizeof(float);
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

void PCMinput::discard(FILE* in, size_t count, const unique_num_array<char>& buffer) const
{	count *= BytesPerSample;
	while (count > buffer.size())
	{	fread2(buffer.get(), buffer.size(), in);
		count -= buffer.size();
	}
	fread2(buffer.get(), count, in);
}
void PCMinput::discard(FILE* in, size_t count)
{	if (!Buffer.size())
		Buffer.reset(BUFFERSIZE);
	discard(in, count, Buffer);
}

void PCMinput::ASCIIdump(FILE* out, const unique_num_array<char>& src)
{	Len = src.size() / BytesPerSample;
	Sp = src.get();
	if (Fmt >= Format::F32)
		while (Len--)
		{	(this->*ReadFunc)();
			fprintf(out, "%g\t%g\n", Ch[0], Ch[1]);
		}
	else
		while (Len--)
		{	(this->*ReadFunc)();
			fprintf(out, "%i\t%i\n", (int)Ch[0], (int)Ch[1]);
		}
}

void PCMinput::read(FILE* in, const unique_num_array<fftw_real>& dstL, const unique_num_array<fftw_real>& dstR, bool add, FILE* asciidump)
{	assert(dstL.size() == dstR.size());
	if (!Buffer.size())
		Buffer.reset(BUFFERSIZE);
	size_t blksize = Buffer.size() / BytesPerSample;
	size_t offset = 0;
	while (dstL.size() - offset > blksize)
	{	fread2(Buffer.begin(), BytesPerSample * blksize, in);
		convert(dstL.slice(offset, blksize), dstR.slice(offset, blksize), Buffer, add);
		if (asciidump)
			ASCIIdump(asciidump, Buffer);
		offset += blksize;
	}
	blksize = dstL.size() - offset;
	fread2(Buffer.begin(), BytesPerSample * blksize, in);
	auto bufslice(Buffer.slice(0, blksize * BytesPerSample));
	convert(dstL.slice(offset, blksize), dstR.slice(offset, blksize), bufslice, add);
	if (asciidump)
		ASCIIdump(asciidump, bufslice);
}
