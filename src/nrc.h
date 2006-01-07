int ludcmp(double* a, int n, int *indx, double* vv);
/* LU decomposition kernel
	a[n][n]	in/out	Matrix data in, LU decomposite out
	n			in			Matrix size
	indx[n]	out		row interchange
	vv[n]		tmp		temporary buffer
	return 				parity
*/

void lubksb(const double *a, int n, const int *indx, double b[]);
/* LU backsubstitution
	a[n][n]	in			LU decomposite
	n			in			Matrix size
	indx[n]	in			row interchange
	b[n]		in/out	RHS vector in, solution vector out
*/

