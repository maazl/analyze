#ifndef MATHX_H_
#define MATHX_H_

#include "nrc.h"
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <array>
#include <complex>
#include <type_traits>


#define M_2PI (2.*M_PI)
#define M_3PI (3.*M_PI)
#ifndef M_PI_2
#define M_PI_2 (M_PI*.5)
#endif
#define M_PI_180 (M_PI/180.)
#define M_180_PI (180./M_PI)


typedef std::complex<double> Complex;


inline double sqr(double v)
{	return v * v;
}
inline unsigned sqr(int v)
{	return v * v;
}
inline unsigned sqr(unsigned v)
{	return v * v;
}
inline uint64_t sqr(int64_t v)
{	return v * v;
}
inline double cbc(double v)
{	return v * v * v;
}

inline double todB(double f)
{	return 20 * std::log10(f);
}
inline double fromdB(double d)
{	return std::pow(10, d / 20);
}

inline double myrand()
{	return std::rand() / (RAND_MAX + 1.);
}

inline double& re(Complex& z)
{	return ((double(&)[2])(z))[0];
}
inline double& im(Complex& z)
{	return ((double(&)[2])(z))[1];
}

template <typename T, size_t N>
using Vector = std::array<T,N>;
template <size_t N>
using VectorD = Vector<double,N>;
template <size_t N>
using VectorC = Vector<Complex,N>;

template <typename T, size_t R, size_t C>
using Matrix = std::array<Vector<T,C>,R>;
template <size_t R, size_t C>
using MatrixD = Matrix<double,R,C>;
template <size_t R, size_t C>
using MatrixC = Matrix<Complex,R,C>;


template <typename T, size_t N>
const T(&operator&(const Vector<T,N>& v))[N]
{	return *(T(*)[N])v.data();
}

template <typename T, size_t N>
Vector<T,N> operator*(const Vector<T,N>& v, T f)
{	Vector<T,N> r;
	T* d = r.begin();
	for (auto& c : v)
		*d++ = c * f;
	return r;
}

template <typename T, size_t N>
VectorD<N> operator/(const VectorD<N>& v, T f)
{	return v * (1./f);
}

template <typename T, size_t N>
T scalar(const Vector<T,N>& s1, const T* s2, size_t inc2)
{	T r = 0;
	const T* s = s1.begin();
	const T* e = s + N;
	do
	{	r += *s++ * *s2;
		s2 += inc2;
	} while (s < e);
	return r;
}

template <typename T, size_t M, size_t N>
Matrix<T,M,N> operator*(const Matrix<T,M,N>& m, T f)
{	Matrix<T,M,N> r;
	const T* s = m.begin()->begin();
	T* d = r.begin()->begin();
	unsigned l = M * N;
	while (l--)
		*d++ = *s++ * f;
	return r;
}
template <typename T, size_t M, size_t N, size_t P>
Matrix<T,M,P> operator*(const Matrix<T,M,N>& m1, const Matrix<T,N,P>& m2)
{	Matrix<T,M,P> r;
	T* d = r.begin();
	for (size_t m = 0; m < M; ++m)
		for (size_t n = 0; n < N; ++n)
			*d++ = scalar(m1[m], &m2[0][n], N);
	return r;
}

template <typename T, size_t N, size_t P>
Vector<T,P> operator*(const Vector<T,N>& v, const Matrix<T,N,P>& m)
{	Vector<T,P> r;
	T* d = r.begin();
	for (size_t n = 0; n < N; ++n)
		*d++ = scalar(v, &m[0][n], N);
	return r;
}



/*template <int N>
double det(const MatrixD<N,N>& m)
{
}*/
template <typename T>
inline T det(const Matrix<T,1,1>& m)
{	return m[0][0];
}
template <typename T>
inline T det(const Matrix<T,2,2>& m)
{	return m[0][0]*m[1][1] - m[0][1]*m[1][0];
}
template <typename T>
inline T det(const Matrix<T,3,3>& m)
{	return m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1]
	     + m[0][1]*m[1][2]*m[2][0] - m[0][1]*m[1][0]*m[2][2]
	     + m[0][2]*m[1][0]*m[2][1] - m[0][2]*m[1][1]*m[2][0];
}

// Not supported for Complex so far.
template <size_t N>
typename std::enable_if<(N > 3), MatrixD<N,N>>::type inverse(const MatrixD<N,N>& m)
{	MatrixD<N,N> s = m;
	MatrixD<N,N> r;
	int indx[N];
	double col[N];
	ludcmp(&s[0][0], N, indx, col);
	for (unsigned j = 0; j < N; j++)
	{	// Find inverse by columns.
		//for (i = 0; i < N; i++)
		//	col[i] = 0.0;
		memset(col, 0, sizeof col);
		col[j] = 1.0;
		lubksb(&s[0][0], N, indx, col);
		for (unsigned i = 0; i < N; i++)
			r[i][j] = col[i];
	}
	return r;
}
template <typename T>
inline Matrix<T,1,1> inverse(const Matrix<T,1,1>& m)
{	return Matrix<T,1,1>{{{1. / det(m)}}};
}
template <typename T>
inline Matrix<T,2,2> inverse(const Matrix<T,2,2>& m)
{	T sca = 1. / det(m);
	return Matrix<T,2,2>{
	{	{ m[1][1]*sca, -m[0][1]*sca }
	,	{ -m[1][0]*sca, m[0][0]*sca }
	} };
}
template <typename T>
inline Matrix<T,3,3> inverse(const Matrix<T,3,3>& m)
{	T sca = 1. / det(m);
	return Matrix<T,3,3>{
	{	{ sca * (m[1][1]*m[2][2] - m[1][2]*m[2][1]), sca * (m[2][1]*m[0][2] - m[2][2]*m[0][1]), sca * (m[0][1]*m[1][2] - m[0][2]*m[1][1]) }
	,	{ sca * (m[1][2]*m[2][0] - m[1][0]*m[2][2]), sca * (m[2][2]*m[0][0] - m[2][0]*m[0][2]), sca * (m[0][2]*m[1][0] - m[0][0]*m[1][2]) }
	,	{ sca * (m[1][0]*m[2][1] - m[1][1]*m[2][0]), sca * (m[2][0]*m[0][1] - m[2][1]*m[0][0]), sca * (m[0][0]*m[1][1] - m[0][1]*m[1][0]) }
	} };
}

#endif // MATHX_H_
