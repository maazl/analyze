#include <stdlib.h>
#include "nrc.h"

template <int M, int N>
class Matrix;

template <int N>
class Vektor
{protected:
	double Data[N];
 public:
 	double operator[](int n) const
 	{	return Data[n];
 	}
 	double& operator[](int n)
 	{	return Data[n];
 	}
/* 	operator Matrix<1,N>&()
 	{	return *(Matrix<1,N>*)*Data;
 	}	
 	operator Matrix<N,1>&()
 	{	return *(Matrix<N,1>*)*Data;
 	}*/	
   const double* GetRaw() const
   {	return Data;
   }
   double* GetRaw()
   {	return Data;
   }
};

template <int N>
Vektor<N> operator*(const Vektor<N>& v, double f)
{	Vektor<N> r;
	double* d = r.GetRaw();
	const double* s = v.GetRaw();
	int i = N;
	do	*d++ = *s++ * f;
	 while (--i);
	return r;
}

template <int N>
Vektor<N> operator/(const Vektor<N>& v, double f)
{	return v * (1./f);
}

template <int M, int N>
class Matrix
{protected:
	double Data[M][N];
 public:
   const Vektor<N>& operator[](int m) const
   {	return *(const Vektor<N>*)Data[m];
   }
   Vektor<N>& operator[](int m)
   {	return *(Vektor<N>*)Data[m];
   }
   const double* GetRaw() const
   {	return *Data;
   }
   double* GetRaw()
   {	return *Data;
   }
};

template <int N>
double scalar(const Vektor<N>& s1, const double* s2, int inc2)
{	double r = 0;
	const double* s = s1.GetRaw();
	const double* e = s + N;
	do
	{	r += *s++ * *s2;
		s2 += inc2;
	} while (s < e);
	return r;
}		 	

template <int M, int N, int P>
Matrix<M,P> operator*(const Matrix<M,N>& m1, const Matrix<N,P>& m2)
{	Matrix<M,P> r;
	double* d = r.GetRaw();
	for (int m = 0; m < M; ++m)
	{	for (int n = 0; n < N; ++n)
		{	*d++ = scalar(m1[m], m2.GetRaw()+n, N);	
	}	}
	return r;			
}

template <int N, int P>
Vektor<P> operator*(const Vektor<N>& v, const Matrix<N,P>& m)
{	Vektor<P> r;
	double* d = r.GetRaw();
	for (int n = 0; n < N; ++n)
	{	*d++ = scalar(v, m.GetRaw()+n, N);	
	}
	return r;			
}

/*template <int N>
double det(const Matrix<N,N>& m)
{	
}*/

double det(const Matrix<1,1>& m)
{	return m[0][0];
}

double det(const Matrix<2,2>& m)
{	return m[0][0]*m[1][1] - m[0][1]*m[1][0];
}

double det(const Matrix<3,3>& m)
{	return m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1]
        + m[0][1]*m[1][2]*m[2][0] - m[0][1]*m[1][0]*m[2][2]
        + m[0][2]*m[1][0]*m[2][1] - m[0][2]*m[1][1]*m[2][0];
}

template <int N>
Matrix<N,N> inverse(const Matrix<N,N>& m)
{	Matrix<N,N> s = m;
	Matrix<N,N> r;
	int indx[N];
	double col[N];
	ludcmp(s.GetRaw(), N, indx, col);
	for (int j = 0; j < N; j++)
	{	// Find inverse by columns.
 		//for (i = 0; i < N; i++)
 		//	col[i] = 0.0;
 		memset(col, 0, sizeof col);
 		col[j] = 1.0;
 		lubksb(s.GetRaw(), N, indx, col);
      for (int i = 0; i < N; i++)
         r[i][j] = col[i];
   }
   return r;
}

Matrix<1,1> inverse(const Matrix<1,1>& m)
{	Matrix<1,1> r;
	r[0][0] = 1. / det(m);
	return r;	
}

Matrix<2,2> inverse(const Matrix<2,2>& m)
{	double sca = 1. / det(m);
	Matrix<2,2> r;
	r[0][0] = m[1][1]*sca;
	r[0][1] = -m[0][1]*sca;
	r[1][0] = -m[1][0]*sca;
	r[1][1] = m[0][0]*sca;
	return r;	
}

Matrix<3,3> inverse(const Matrix<3,3>& m)
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

template <int O>
class PCA
{
 protected:
 	double Data[O*(O+1)/2];
 	double Wsum;
 	int N;
 public:
 	PCA()
 	{	memset(Data, 0, sizeof Data);
 		Wsum = 0;
 		N = 0;
 	}
 	void Store(const double(&v)[O], double w = 1);
	Vektor<O-1> Result() const;
};

template <int O>
void PCA<O>::Store(const double(&v)[O], double w)
{	const double* e = v + O;
	const double* v1 = v;
	double* d = Data;
	do
	{	const double* v2 = v1;
		do	*d++ += w * *v1 * *v2;				
		 while (++v2 < e);
	} while (++v1 < e);
	Wsum += w;
	++N;
}

template <int O>
Vektor<O-1> PCA<O>::Result() const
{	// design matrix
	Matrix<O-1,O-1> dm;
	const double* s = Data + O;
	for (int m = 0; m < O-1; ++m)
	{	dm[m][m] = *s++;
		for (int n = m+1; n < O-1; ++n)
			dm[m][n] = dm[n][m] = *s++;
	}
	// inverse
 	dm = inverse(dm);
 	return *(const Vektor<O-1>*)&Data[1] * dm;
}

