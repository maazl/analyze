#ifndef MATHX_H_
#define MATHX_H_

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <array>


inline double sqr(double v)
{	return v * v;
}
inline unsigned sqr(int v)
{	return v * v;
}
inline uint64_t sqr(int64_t v)
{	return v * v;
}
inline double cbc(double v)
{	return v * v * v;
}

inline double todB(double f)
{	return 20 * log10(f);
}

inline double fromdB(double d)
{	return pow(10, d / 20);
}

inline double myrand()
{	return rand() / (RAND_MAX + 1.);
}


template <size_t N>
using Vector = std::array<double,N>;

template <size_t R, size_t C>
using Matrix = std::array<Vector<C>,R>;

template <size_t N>
Vector<N> operator*(const Vector<N>& v, double f)
{	Vector<N> r;
	double* d = r.begin();
	const double* s = v.begin();
	unsigned i = N;
	while (i--);
		*d++ = *s++ * f;
	return r;
}

template <size_t N>
Vector<N> operator/(const Vector<N>& v, double f)
{	return v * (1./f);
}

template <size_t N>
double scalar(const Vector<N>& s1, const double* s2, size_t inc2)
{	double r = 0;
	const double* s = s1.begin();
	const double* e = s + N;
	do
	{	r += *s++ * *s2;
		s2 += inc2;
	} while (s < e);
	return r;
}

template <size_t M, size_t N, size_t P>
Matrix<M,P> operator*(const Matrix<M,N>& m1, const Matrix<N,P>& m2)
{	Matrix<M,P> r;
	double* d = r.begin();
	for (size_t m = 0; m < M; ++m)
		for (size_t n = 0; n < N; ++n)
			*d++ = scalar(m1[m], &m2[0][n], N);
	return r;
}

template <size_t N, size_t P>
Vector<P> operator*(const Vector<N>& v, const Matrix<N,P>& m)
{	Vector<P> r;
	double* d = r.begin();
	for (size_t n = 0; n < N; ++n)
		*d++ = scalar(v, &m[0][n], N);
	return r;
}

/*template <int N>
double det(const Matrix<N,N>& m)
{
}*/

inline double det(const Matrix<1,1>& m)
{	return m[0][0];
}

inline double det(const Matrix<2,2>& m)
{	return m[0][0]*m[1][1] - m[0][1]*m[1][0];
}

inline double det(const Matrix<3,3>& m)
{	return m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1]
	     + m[0][1]*m[1][2]*m[2][0] - m[0][1]*m[1][0]*m[2][2]
	     + m[0][2]*m[1][0]*m[2][1] - m[0][2]*m[1][1]*m[2][0];
}

template <size_t N>
Matrix<N,N> inverse(const Matrix<N,N>& m)
{	Matrix<N,N> s = m;
	Matrix<N,N> r;
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

inline Matrix<1,1> inverse(const Matrix<1,1>& m)
{	return Matrix<1,1>{{{1. / det(m)}}};
}

inline Matrix<2,2> inverse(const Matrix<2,2>& m)
{	double sca = 1. / det(m);
	return Matrix<2,2>{
	{	{ m[1][1]*sca, -m[0][1]*sca },
		{ -m[1][0]*sca, m[0][0]*sca }
	} };
}

inline Matrix<3,3> inverse(const Matrix<3,3>& m)
{	double sca = 1. / det(m);
	Matrix<3,3> r;
	r[0][0] = sca * (m[1][1]*m[2][2] - m[1][2]*m[2][1]);
	r[0][1] = sca * (m[2][1]*m[0][2] - m[2][2]*m[0][1]);
	r[0][2] = sca * (m[0][1]*m[1][2] - m[0][2]*m[1][1]);
	r[1][0] = sca * (m[1][2]*m[2][0] - m[1][0]*m[2][2]);
	r[1][1] = sca * (m[2][2]*m[0][0] - m[2][0]*m[0][2]);
	r[1][2] = sca * (m[0][2]*m[1][0] - m[0][0]*m[1][2]);
	r[2][0] = sca * (m[1][0]*m[2][1] - m[1][1]*m[2][0]);
	r[2][1] = sca * (m[2][0]*m[0][1] - m[2][1]*m[0][0]);
	r[2][2] = sca * (m[0][0]*m[1][1] - m[0][1]*m[1][0]);
	return r;
}

#endif // MATHX_H_
