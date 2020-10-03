#ifndef SCOPED_ARRAY_H_
#define SCOPED_ARRAY_H_

#include <cstddef>
#include <memory>
#include <cstring>
#include <cassert>

/*template<class T> class std::complex;
template<class T> struct is_complex : std::false_type {};
template<class T> struct is_complex<std::complex<T>> : std::true_type {};*/

/// Array of T with ownership. Like \see std::unique_ptr<T[]> but with size tracking.
/// @tparam T Element type.
template <typename T>
class unique_array : public std::unique_ptr<T[],void (*)(void*)>
{	size_t Size;
 private:
	typedef std::unique_ptr<T[],void (*)(void*)> base;
 protected:
	unique_array(T* ptr, size_t size, void (*deleter)(void*)) noexcept : base(ptr, deleter), Size(size) {}
	void reset(T* ptr, size_t size) noexcept { base::reset(ptr); Size = size; }
 public:
	constexpr unique_array() noexcept : base(nullptr, operator delete[]), Size(0) {}
	explicit unique_array(size_t size) : base(new T[size], operator delete[]), Size(size) {}
	unique_array(unique_array<T>&& r) : base(move(r)), Size(r.Size) { r.Size = 0; }
	void reset(size_t size = 0) { base::reset(size ? new T[size] : nullptr); Size = size; }
	void swap(unique_array<T>&& other) noexcept { base::swap(other); std::swap(Size, other.Size); }
	void assign(const unique_array<T>& r) const { assert(this->size() == r.size()); memmove(this->begin(), r.begin(), this->size() * sizeof(T)); }
	const unique_array<T>& operator =(const unique_array<T>& r) const { assign(r); return *this; }
	size_t size() const noexcept { return Size; }
	T* begin() const noexcept { return base::get(); }
	T* end() const noexcept { return begin() + Size; }
	T& operator[](size_t i) const { assert(i < Size); return base::operator[](i); }
	unique_array<T> slice(size_t start, size_t count) const noexcept
	{	assert(start + count <= Size); return unique_array<T>(begin() + start, count, [](void*){}); }
};

/// Array of numeric T with ownership and mathematical vector operations.
/// @tparam T Numeric element type.
template <typename T>
class unique_num_array : public unique_array<T>
{	// Neither is_arithmetic nor is_standard_layout holds for std::complex
	//static_assert(std::is_arithmetic<T>::value, "T must be arithmetic");
 public:
	using unique_array<T>::unique_array;
	unique_num_array(unique_num_array<T>&& r) : unique_array<T>(move(r)) {}
	unique_num_array() {}
	unique_num_array<T> slice(size_t start, size_t count) const noexcept
	{	assert(start + count <= this->size()); return unique_num_array<T>(this->begin() + start, count, [](void*){}); }
 public: // math operations, require some numeric type T
	void clear() const { std::memset(this->begin(), 0, this->size() * sizeof(T)); }
	const unique_num_array<T>& operator =(const unique_num_array<T>& r) const { this->assign(r); return *this; }
	const unique_num_array<T>& operator +=(const unique_num_array<T>& r) const;
	const unique_num_array<T>& operator -=(const unique_num_array<T>& r) const;
	const unique_num_array<T>& operator *=(T r) const;
	const unique_num_array<T>& operator *=(const unique_num_array<T>& r) const;
	const unique_num_array<T>& operator /=(T r) const { return *this *= 1/r; }
	const unique_num_array<T>& operator /=(const unique_num_array<T>& r) const;
	//std::enable_if<is_complex<T>, const unique_num_array<T>&> conj() const;
};

template <typename T>
const unique_num_array<T>& unique_num_array<T>::operator +=(const unique_num_array<T>& r) const
{	size_t len = this->size();
	assert(len == r.size());
	T* dp = this->begin();
	const T* sp = r.begin();
	while (len--)
		*dp++ += *sp++;
	return *this;
}

template <typename T>
const unique_num_array<T>& unique_num_array<T>::operator -=(const unique_num_array<T>& r) const
{	size_t len = this->size();
	assert(len == r.size());
	T* dp = this->begin();
	const T* sp = r.begin();
	while (len--)
		*dp++ -= *sp++;
	return *this;
}

template <typename T>
const unique_num_array<T>& unique_num_array<T>::operator *=(T r) const
{	size_t len = this->size();
	T* dp = this->begin();
	while (len--)
		*dp++ *= r;
	return *this;
}
template <typename T>
const unique_num_array<T>& unique_num_array<T>::operator *=(const unique_num_array<T>& r) const
{	size_t len = this->size();
	assert(len == r.size());
	T* dp = this->begin();
	const T* sp = r.begin();
	while (len--)
		*dp++ *= *sp++;
	return *this;
}

template <typename T>
const unique_num_array<T>& unique_num_array<T>::operator /=(const unique_num_array<T>& r) const
{	size_t len = this->size();
	assert(len == r.size());
	T* dp = this->begin();
	const T* sp = r.begin();
	while (len--)
		*dp++ /= *sp++;
	return *this;
}

/*template <typename T>
std::enable_if<is_complex<T>, const unique_num_array<T>&> unique_num_array<T>::conj() const
{	T* dp = this->begin();
	for (T* dp = this->begin(), ep = dp + this->size(); dp != ep; ++dp)
		*dp = std::conj(*dp);
	return *this;
}*/

extern "C"
{	void *fftwf_malloc(size_t n);
	void fftwf_free(void *p);
}
template <typename T>
class unique_fftw_arr : public unique_num_array<T>
{private:
	typedef unique_num_array<T> base;
 public:
	unique_fftw_arr() noexcept : base(NULL, 0, fftwf_free) {}
	explicit unique_fftw_arr(size_t size) : base((T*)fftwf_malloc(size * sizeof(T)), size, fftwf_free) {}
	void reset(size_t size = 0) { base::reset(size ? (T*)fftwf_malloc(size * sizeof(T)) : NULL, size); }
};


#endif // SCOPED_ARRAY_H_
