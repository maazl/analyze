#include "moment.h"

void Skewness::push(double val, double weight)
{	Sum3 += cb(val) * weight; Variance::push(val, weight);
}

double Skewness::skewness() const
{	return (Sum3 - (3 * Sum2 * Sum - 2 * cb(Sum) / Count) / Count) / Count / pow(varianceMLE(), 1.5);
}

double Kurtosis::kurtosis() const
{	return (Sum4 - 4 * Sum3 * Sum / Count + 6 * Sum2 * sqr(Sum / Count) - 3 * sqr(sqr(Sum)) / cb(Count)) / Count / sqr(varianceMLE());
}
double Kurtosis::kurtosisExcessBC() const
{	return (kurtosis() * (Count + 1) - 3 * (Count - 1)) * (Count - 1) / ((Count - 2) * (Count - 3));
}
