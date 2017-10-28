#include "utils.h"

#include <stdarg.h>
#include <stdlib.h>
#include <inttypes.h>


bool termrq = false;

void die(int rc, const char* msg, ...)
{
	termrq = true;
	va_list va;
	va_start(va, msg);
	vfprintf(stderr, msg, va);
	va_end(va);
	fputc('\n', stderr);
	exit(rc);
}

FILE* binmode(FILE* stream)
{	// TODO: Windows & Co ...
	return stream;
}

void wavheader(FILE* fo, size_t nsamp, size_t sfreq)
{
	uint32_t wavhdr[11] =
	{	0x46464952, 0xffffffff, 0x45564157, 0x20746D66,
		0x00000010, 0x00020001, 0x0000BB80, 0x0002ee00,
		0x00100004, 0x61746164, 0xffffffff };

	wavhdr[6] = sfreq;
	wavhdr[7] = sfreq * 4;

	wavhdr[10] = nsamp * sizeof(short);
	wavhdr[1] = wavhdr[10] + 44;

	fwrite(wavhdr, sizeof wavhdr, 1, fo);
}

void fread2(void* data, size_t size, size_t count, FILE* stream)
{	while (count)
	{	size_t read = fread(data, size, count, stream);
		if (read == 0)
			die(27, "Failed to read data from input.");
		(char*&)data += size * read;
		count -= read;
	}
}
