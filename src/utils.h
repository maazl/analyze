#ifndef UTILS_H_
#define UTILS_H_

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <array>

#ifdef __GNUC__
#define PRINTFATTR(i) __attribute__((format(printf, i, i+1)))
#define SCANFATTR(i) __attribute__((format(scanf, i, i+1)))
#else
#define PRINTFATTR(i)
#define SCANFATTR(i)
#endif


// Termination flag, set by die(...).
extern bool termrq;

/** Terminate program with exit code and error message
 * and set termrq to notify other threads (if any).
 * @param rc exit code
 * @param msg error message format string
 * @param ... error message arguments
 */
void die(int rc, const char* msg, ...) PRINTFATTR(2);

/** Switch stream to binary mode
 * @param stream file stream
 * @return modified file stream - might be the same than \a stream or not.
 */
std::FILE* binmode(std::FILE* stream);

/// Checked file open. Like fopen but dies on error.
/// @param file Name of the file, "-" = stdin/stdout
/// @param mode Open mode
/// @return File pointer, never NULL.
/// @exception Message Failed to open: error message
std::FILE* checkedopen(const char* file, const char* mode);

/// RAII version of FILE*
class FILEguard
{	std::FILE* File;
 public:
	FILEguard(std::FILE* file = nullptr) : File(file) {}
	FILEguard(const FILEguard&) = delete;
	FILEguard(FILEguard&& r) : File(r.File) { r.File = nullptr; }
	FILEguard(const char* file, const char* mode) : File(checkedopen(file, mode)) {}
	~FILEguard() { if (File) std::fclose(File); }
	FILEguard& operator=(const FILEguard&) = delete;
	FILEguard& operator=(FILE* file) { this->~FILEguard(); File = file; return *this; }
	operator std::FILE*() const { return File; }
};

/** Write RIFF WAV header.
 * @param fo output stream
 * @param nsamp number of samples to write
 * @param sfreq sampling frequency
 */
void wavheader(std::FILE* fo, size_t nsamp, size_t sfreq);

/** Really read count items or die.
 * @param data target buffer
 * @param count number of bytes
 * @param stream source stream
 */
void fread2(void* data, size_t count, std::FILE* stream);

/** Really write count items or die.
 * @param data target buffer
 * @param count number of bytes
 * @param stream source stream
 */
void fwrite2(const void* data, size_t count, std::FILE* stream);


extern const int16_t endian_detect_;
static inline bool is_big_endian()
{	return (bool)*(char*)&endian_detect_;
}

inline uint16_t bswap(uint16_t v)
{	//return _srotl(v, 8);
	return (uint16_t)v >> 8 | v << 8;
}

template<size_t N>
inline void cswap(char* cp)
{	std::swap(cp[0], cp[N-1]);
	cswap<N-2>(cp+1);
}
template<>
inline void cswap<1>(char*)
{}
template<>
inline void cswap<0>(char*)
{}
template<typename T>
inline void bswap(T* dp)
{	cswap<sizeof(T)>((char*)dp);
}

/** Execute shell command or die.
 * @param cmd Command to execute, if \c nullptr then this is a no-op.
 */
void execute(const char* cmd);


/// Get rid of GCC's unused result warnings.
template<typename T>
inline void ignore_result(T)
{}

/// Reference helper to allow array of references
template <typename T>
struct reference
{	T& Ref;
	constexpr operator T&() const { return Ref; }
};

/// learn integer_sequence if < C++14
#if __cplusplus < 201402L
template<std::size_t... Is> struct integer_sequence{};
template<std::size_t N, std::size_t... Is>
struct make_index_sequence : make_index_sequence<N-1, N-1, Is...>{};
template<std::size_t... Is>
struct make_index_sequence<0, Is...> : integer_sequence<Is...>{};
#else
template<std::size_t... Is>
using integer_sequence = std::integer_sequence<Is...>;
template<std::size_t N>
using make_index_sequence = std::make_integer_sequence<std::size_t, N>;
#endif

/// constexpr string, basically std::array<const char,N>,
/// but constructible from string literal of different length.
template <std::size_t N>
class cstring : public std::array<const char,N>
{private:
	template<std::size_t M, std::size_t... I>
	constexpr cstring(const char (&v)[M], integer_sequence<I...>) : std::array<const char,N>{{ v[I]... }} {}
 public:
	/// Create from string literal
	template <std::size_t M>
	constexpr cstring(const char (&v)[M]) : cstring(v, make_index_sequence<M-1>{}) {}
};

#endif // UTILS_H_
