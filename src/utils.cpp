#include "utils.h"

#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <cinttypes>
#include <cerrno>

using namespace std;

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

FILE* checkedopen(const char* file, const char* mode)
{	FILE* ret;
	if (strcmp(file, "-") == 0)
	{	if (strchr(mode, 'r'))
			ret = stdin;
		else if (strchr(mode, 'w') || strrchr(mode, 'a'))
			ret = stdout;
		else
			die(28, "invalid open mode.");
		if (strchr(mode, 'b'))
			ret = binmode(ret);
	} else
	{	ret = fopen(file, mode);
		if (ret == NULL)
			die(21, "Failed to open file '%s' for %s: %s", file, strchr(mode, 'r') ? "reading" : "writing", strerror(errno));
	}
	return ret;
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

void fread2(void* data, size_t count, FILE* stream)
{	while (count)
	{	size_t read = fread(data, 1, count, stream);
		if (read == 0)
			die(27, "Failed to read %zu bytes: %s", count, ferror(stream) ? strerror(errno) : "end of file");
		count -= read;
		(char*&)data += read;
	}
}

void fwrite2(const void* buffer, size_t count, FILE* stream)
{	while (count)
	{	size_t wrote = fwrite(buffer, 1, count, stream);
		if (wrote == 0)
			die(27, "Failed to write %zu bytes: %s", count, ferror(stream) ? strerror(errno) : "no space left");
		count -= wrote;
		(const char*&)buffer += wrote;
	}
}

const int16_t endian_detect_ = 0x0100;

void execute(const char* cmd)
{	if (cmd)
	{	int rc = system(cmd);
		if (rc)
			die(25, "Failed to execute external command '%.30s', RC = %i.", cmd, rc);
	}
}
