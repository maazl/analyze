#include "pcmio.h"
#include "mathx.h"
#include "utils.h"

#include <float.h>
#include <assert.h>

using namespace std;

#define BUFFERSIZE 65536

PCMoutput::PCMoutput(Format fmt, fftw_real gain, bool symmetric)
:	PCMIO(fmt)
,	Gain(fmt >= Format::F32 ? gain : gain * (getFSR(fmt) - 1))
,	Symmetric(symmetric)
{	fprintf(stderr, "G: %g\t%g\t%g\t%i\n", Gain, gain, getFSR(fmt), (int)fmt);
}

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
				*dp++ = Symmetric ? -s : s;
		}	}
		break;
	 case Format::F32_SWAP:
		{	float* dp = (float*)dst.begin();
			while (sp != spe)
			{	dp[0] = *sp++ * Gain;
				if (Symmetric)
				{	dp[1] = -dp[0];
					bswap(dp);
					bswap(dp + 1);
				} else
				{	bswap(dp);
					dp[1] = dp[0];
				}
				dp += 2;
		}	}
		break;
	 case Format::I32:
		{	int32_t* dp = (int32_t*)dst.begin();
			while (sp != spe)
			{	dp[0] = (int32_t)round(*sp++ * Gain);
				dp[1] = Symmetric ? -dp[0] : dp[0];
				dp += 2;
		}	}
		break;
	 case Format::I32_SWAP:
		{	int32_t* dp = (int32_t*)dst.begin();
			while (sp != spe)
			{	dp[0] = round(*sp++ * Gain);
				if (Symmetric)
				{	dp[1] = -dp[0];
					bswap(dp);
					bswap(dp + 1);
				} else
				{	bswap(dp);
					dp[1] = dp[0];
				}
				dp += 2;
		}	}
		break;
	 case Format::I24:
		{	char* dp = (char*)dst.begin();
			endian_conv e;
			char* cp = e.c + is_big_endian();
			while (sp != spe)
			{	e.i = (int32_t)round(*sp++ * Gain);
				dp[0] = cp[0];
				dp[1] = cp[1];
				dp[2] = cp[2];
				if (Symmetric)
					e.i = -e.i;
				dp[3] = cp[0];
				dp[4] = cp[1];
				dp[5] = cp[2];
				dp += 6;
		}	}
		break;
	 case Format::I24_SWAP:
		{	char* dp = (char*)dst.begin();
			endian_conv e;
			char* cp = e.c + is_big_endian();
			while (sp != spe)
			{	e.i = (int32_t)round(*sp++ * Gain);
				dp[0] = cp[2];
				dp[1] = cp[1];
				dp[2] = cp[0];
				if (Symmetric)
					e.i = -e.i;
				dp[3] = cp[2];
				dp[4] = cp[1];
				dp[5] = cp[0];
				dp += 6;
		}	}
		break;
	 case Format::I16:
		{	int16_t* dp = (int16_t*)dst.begin();
			while (sp != spe)
			{	dp[0] = (int16_t)floor(*sp++ * Gain + myrand());
				dp[1] = Symmetric ? -dp[0] : dp[0];
				dp += 2;
		}	}
		break;
	 case Format::I16_SWAP:
		{	int16_t* dp = (int16_t*)dst.begin();
			if (Symmetric)
				while (sp != spe)
				{	int16_t s = (int16_t)floor(*sp++ * Gain + myrand());
					*dp++ = bswap(s);
					*dp++ = bswap(-s);
				}
			else
				while (sp != spe)
				{	dp[0] = dp[1] = bswap((int16_t)floor(*sp++ * Gain + myrand()));
					dp += 2;
		}		}
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
	 case Format::F32_SWAP:
		{	float* dp = (float*)dst.begin();
			while (sp1 != spe)
			{	*dp = *sp1++ * Gain;
				bswap(dp++);
				*dp = *sp2++ * Gain;
				bswap(dp++);
				dp += 2;
		}	}
		break;
	 case Format::I32:
		{	int32_t* dp = (int32_t*)dst.begin();
			while (sp1 != spe)
			{	*dp++ = (int32_t)round(*sp1++ * Gain);
				*dp++ = (int32_t)round(*sp2++ * Gain);
		}	}
		break;
	 case Format::I32_SWAP:
		{	int32_t* dp = (int32_t*)dst.begin();
			while (sp1 != spe)
			{	*dp = (int32_t)round(*sp1++ * Gain);
				bswap(dp++);
				*dp = (int32_t)round(*sp2++ * Gain);
				bswap(dp++);
				dp += 2;
		}	}
		break;
	 case Format::I24:
		{	char* dp = dst.begin();
			endian_conv e;
			char* cp = e.c + is_big_endian();
			while (sp1 != spe)
			{	e.i = (int32_t)round(*sp1++ * Gain);
				dp[0] = cp[0];
				dp[1] = cp[1];
				dp[2] = cp[2];
				e.i = (int32_t)round(*sp2++ * Gain);
				dp[3] = cp[0];
				dp[4] = cp[1];
				dp[5] = cp[2];
				dp += 6;
		}	}
		break;
	 case Format::I24_SWAP:
		{	char* dp = dst.begin();
			endian_conv e;
			char* cp = e.c + is_big_endian();
			while (sp1 != spe)
			{	e.i = (int32_t)round(*sp1++ * Gain);
				dp[0] = cp[2];
				dp[1] = cp[1];
				dp[2] = cp[0];
				e.i = (int32_t)round(*sp2++ * Gain);
				dp[3] = cp[2];
				dp[4] = cp[1];
				dp[5] = cp[0];
				dp += 6;
		}	}
		break;
	 case Format::I16:
		{	int16_t* dp = (int16_t*)dst.begin();
			while (sp1 != spe)
			{	*dp++ = (int16_t)floor(*sp1++ * Gain + myrand());
				*dp++ = (int16_t)floor(*sp2++ * Gain + myrand());
		}	}
		break;
	 case Format::I16_SWAP:
		{	int16_t* dp = (int16_t*)dst.begin();
			while (sp1 != spe)
			{	*dp++ = bswap((int16_t)floor(*sp1++ * Gain + myrand()));
				*dp++ = bswap((int16_t)floor(*sp2++ * Gain + myrand()));
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
