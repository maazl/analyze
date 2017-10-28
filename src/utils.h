#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>

// Termination flag, set by die(...).
extern bool termrq;

/** Terminate program with exit code and error message
 * and set termrq to notify other threads (if any).
 * @param rc exit code
 * @param msg error message format string
 * @param ... error message arguments
 */
void die(int rc, const char* msg, ...);

/** Switch stream to binary mode
 * @param stream file stream
 * @return modified file stream - might be the same than \a stream or not.
 */
FILE* binmode(FILE* stream);

/** Write RIFF WAV header.
 * @param fo output stream
 * @param nsamp number of samples to write
 * @param sfreq sampling frequency
 */
void wavheader(FILE* fo, size_t nsamp, size_t sfreq);

/** Really read up to count items (or die).
 * @param data target buffer
 * @param size size of items
 * @param count number of items
 * @param stream source stream
 */
void fread2(void* data, size_t size, size_t count, FILE* stream);

#endif // UTILS_H_
