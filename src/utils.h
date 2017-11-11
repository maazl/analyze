#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <memory>


template <typename T>
class scoped_array : public std::unique_ptr<T[],void (*)(void*)>
{	//static_assert(std::is_pod<T>::value, "T must be POD");
	size_t Size;
 protected:
	typedef std::unique_ptr<T[],void (*)(void*)> base;
	scoped_array(T* ptr, size_t size, void (*deleter)(void*)) noexcept : base(ptr, deleter), Size(size) {}
	void reset(T* ptr, size_t size) noexcept { base::reset(ptr); Size = size; }
 public:
	constexpr scoped_array() noexcept : base(NULL, operator delete[]), Size(0) {}
	explicit scoped_array(size_t size) : base(new T[size], operator delete[]), Size(size) {}
	void reset(size_t size = 0) { base::reset(size ? new T[size] : NULL); Size = size; }
	void swap(scoped_array<T>& other) noexcept { base::swap(other); swap(Size, other.Size); }
	void clear() const { memset(base::get(), 0, Size * sizeof(T)); }
	void copyfrom(const scoped_array<T>& other) const { assert(Size == other.Size); memcpy(base::get(), other.get(), Size * sizeof(T)); }
	size_t size() const noexcept { return Size; }
	T* begin() const noexcept { return base::get(); }
	T* end() const noexcept { return base::get() + Size; }
	T& operator[](size_t i) const { assert(i < Size); return base::operator[](i); }
	scoped_array<T> slice(size_t start, size_t count) const noexcept
	{	return scoped_array<T>(begin() + start, count, [](void*){}); }
};

extern "C"
{	void *fftw_malloc(size_t n);
	void fftw_free(void *p);
}
template <typename T>
class scoped_fftw_arr : public scoped_array<T>
{protected:
	typedef scoped_array<T> base;
 public:
	scoped_fftw_arr() noexcept : base(NULL, 0, fftw_free) {}
	explicit scoped_fftw_arr(size_t size) : base((T*)fftw_malloc(size * sizeof(T)), size, fftw_free) {}
	void reset(size_t size = 0) { base::reset(size ? (T*)fftw_malloc(size * sizeof(T)) : NULL, size); }
};


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

/// Checked file open. Like fopen but throws on error.
/// @param file Name of the file
/// @param mode Open mode
/// @return File pointer, never NULL.
/// @exception Message Failed to open: error message
FILE* checkedopen(const char* file, const char* mode);

/// RAII version of FILE*
class FILEguard
{	FILE* File;
 public:
	FILEguard(FILE* file) : File(file) {}
	FILEguard(const FILEguard&) = delete;
	FILEguard(FILEguard&& r) : File(r.File) { r.File = NULL; }
	FILEguard(const char* file, const char* mode) : File(checkedopen(file, mode)) {}
	~FILEguard() { if (File) fclose(File); }
	FILEguard& operator=(const FILEguard&) = delete;
	FILEguard& operator=(FILE* file) { this->~FILEguard(); File = file; return *this; }
	operator FILE*() { return File; }
};

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

static inline uint16_t bswap(uint16_t v)
{	//return _srotl(v, 8);
	return (uint16_t)v >> 8 | v << 8;
}

static inline double sqr(double v)
{	return v * v;
}
static inline int64_t sqr(int64_t v)
{	return v * v;
}
static inline double cb(double v)
{	return v * v * v;
}

static inline double todB(double f)
{	return 20 * log10(f);
}

static inline double fromdB(double d)
{	return pow(10, d / 20);
}


#endif // UTILS_H_
