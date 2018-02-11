#include "pcmio.h"
#include "mathx.h"
#include "utils.h"

#include <float.h>
#include <assert.h>

using namespace std;

#define BUFFERSIZE 65536


PCMoutput::PCMoutput(Format fmt, fftw_real gain)
:	PCMIO(fmt)
,	Gain(fmt == Format::F32 ? gain : 32767 * gain)
{}

void PCMoutput::convert(const unique_num_array<fftw_real>& src, const unique_num_array<char>& dst)
{	assert(BytesPerSample * src.size() == dst.size());
	const fftw_real* sp = src.begin();
	const fftw_real* const spe = sp + src.size();
	switch (Fmt)
	{case Format::F32:
		{	float* dp = (float*)dst.begin();
			while (sp != spe)
			{	float s = *sp++ * Gain;
				*dp++ = s;
				*dp++ = -s;
		}	}
		break;
	 case Format::I16:
		{	short* dp = (short*)dst.begin();
			while (sp != spe)
			{	short s = (short)floor(*sp++ * Gain + myrand());
				*dp++ = s;
				*dp++ = -s;
		}	}
		break;
	 case Format::I16_SWAP:
		{	short* dp = (short*)dst.begin();
			while (sp != spe)
			{	short s = bswap((short)floor(*sp++ * Gain + myrand()));
				*dp++ = s;
				*dp++ = -s;
		}	}
		break;
	}
}
void PCMoutput::convert(const unique_num_array<fftw_real>& src1, const unique_num_array<fftw_real>& src2, const unique_num_array<char>& dst)
{	assert(BytesPerSample * src1.size() == dst.size() && src1.size() == src2.size());
	const fftw_real* sp1 = src1.begin();
	const fftw_real* const spe = sp1 + src1.size();
	const fftw_real* sp2 = src2.begin();
	switch (Fmt)
	{case Format::F32:
		{	float* dp = (float*)dst.begin();
			while (sp1 != spe)
			{	*dp++ = *sp1++ * Gain;
				*dp++ = *sp2++ * Gain;
		}	}
		break;
	 case Format::I16:
		{	short* dp = (short*)dst.begin();
			while (sp1 != spe)
			{	*dp++ = (short)floor(*sp1++ * Gain + myrand());
				*dp++ = (short)floor(*sp2++ * Gain + myrand());
		}	}
		break;
	 case Format::I16_SWAP:
		{	short* dp = (short*)dst.begin();
			while (sp1 != spe)
			{	*dp++ = bswap((short)floor(*sp1++ * Gain + myrand()));
				*dp++ = bswap((short)floor(*sp2++ * Gain + myrand()));
		}	}
		break;
	}
}

void PCMoutput::write(FILE* out, const unique_num_array<fftw_real>& src)
{	if (!Buffer.size())
		Buffer.reset(BUFFERSIZE);
	size_t blksize = Buffer.size() / BytesPerSample;
	size_t offset = 0;
	while (src.size() - offset > blksize)
	{	convert(src.slice(offset, blksize), Buffer);
		fwrite2(Buffer.begin(), BytesPerSample * blksize, out);
		offset += blksize;
	}
	blksize = src.size() - offset;
	convert(src.slice(offset, blksize), Buffer.slice(0, blksize * BytesPerSample));
	fwrite2(Buffer.begin(), BytesPerSample * blksize, out);
}
void PCMoutput::write(FILE* out, const unique_num_array<fftw_real>& srcL, const unique_num_array<fftw_real>& srcR)
{	assert(srcL.size() == srcR.size());
	if (!Buffer.size())
		Buffer.reset(BUFFERSIZE);
	size_t blksize = Buffer.size() / BytesPerSample;
	size_t offset = 0;
	while (srcL.size() - offset > blksize)
	{	convert(srcL.slice(offset, blksize), srcR.slice(offset, blksize), Buffer);
		fwrite2(Buffer.begin(), BytesPerSample * blksize, out);
		offset += blksize;
	}
	blksize = srcL.size() - offset;
	convert(srcL.slice(offset, blksize), srcR.slice(offset, blksize), Buffer.slice(0, blksize * BytesPerSample));
	fwrite2(Buffer.begin(), BytesPerSample * blksize, out);
}

void PCMoutput::zero(FILE* out, size_t count)
{	if (!Buffer.size())
		Buffer.reset(BUFFERSIZE);
	Buffer.clear();
	count *= BytesPerSample;
	while (count > Buffer.size())
	{	fwrite2(Buffer.begin(), Buffer.size(), out);
		count -= Buffer.size();
	}
	fwrite2(Buffer.begin(), count, out);
}

void PCMoutput::WAVheader(FILE* fo, size_t nsamp, size_t sfreq)
{
	uint32_t wavhdr[11] =
	{	0x46464952, 0xffffffff, 0x45564157, 0x20746D66,
		0x00000010, 0x00020001, 0x0000BB80, 0x0002ee00,
		0x00100004, 0x61746164, 0xffffffff };

	wavhdr[6] = sfreq;
	wavhdr[7] = sfreq * BytesPerSample;

	wavhdr[10] = nsamp ? nsamp * BytesPerSample : (0x7fffffff-sizeof wavhdr) & -BytesPerSample;
	wavhdr[1] = wavhdr[10] + 36;

	fwrite2(wavhdr, sizeof wavhdr, fo);
}
