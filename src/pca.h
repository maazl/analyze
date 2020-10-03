#include "nrc.h"
#include "mathx.h"
#include <stdlib.h>
#include <string.h>


/// Class to do Principal Component Analysis.
/// @tparam O Number of components + 1.
template <size_t O>
class PCA
{protected:
	double Data[O*(O+1)/2];
	double Wsum;
	size_t N;
 public:
	void reset()
	{	memset(Data, 0, sizeof Data);
		Wsum = 0;
		N = 0;
	}
	PCA()
	{	reset();
	}
	/// Add a single data tuple to the PCA buffer. You can add as many data points as you want,
	/// but you will need at least O data points for the result to become defined.
	/// @param v Vector with data value (the first element) and the corresponding values of the components.
	/// @param w Weight of this value. 1 by default.
	void Store(const double(&v)[O], double w = 1);
	/// Calculate the result.
	/// @return Vector with factors for each component. Using the component function with this factors
	/// will reproduce the result as accurate as possible.
	VectorD<O-1> Result() const;
};

template <size_t O>
void PCA<O>::Store(const double(&v)[O], double w)
{	const double* e = v + O;
	const double* v1 = v;
	double* d = Data;
	do
	{	const double* v2 = v1;
		do
			*d++ += w * *v1 * *v2;
		while (++v2 < e);
	} while (++v1 < e);
	Wsum += w;
	++N;
}

template <size_t O>
VectorD<O-1> PCA<O>::Result() const
{	// design matrix
	MatrixD<O-1,O-1> dm;
	const double* s = Data + O;
	for (unsigned m = 0; m < O-1; ++m)
	{	dm[m][m] = *s++;
		for (unsigned n = m+1; n < O-1; ++n)
			dm[m][n] = dm[n][m] = *s++;
	}
	// inverse
 	dm = inverse(dm);
 	return *(const VectorD<O-1>*)&Data[1] * dm;
}

