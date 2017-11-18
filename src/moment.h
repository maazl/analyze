#ifndef MOMENT_H_
#define MOMENT_H_

#include "math.h"
#include "mathx.h"

class Mean
{protected:
	double Count;
	double Sum;
 public:
	Mean() : Count(0), Sum(0) {}
	void push(double val, double weight = 1) { Count += weight; Sum += val * weight; }
	double count() const { return Count; }
	double sum() const { return Sum; }
  double mean() const { return Sum / Count; }
};

class Variance : public Mean
{protected:
	double Sum2;
	//double varsum() const { return (Sum2 - sqr(Sum) / Count); }
 public:
	Variance() : Sum2(0) {}
	void push(double val, double weight = 1) { Sum2 += sqr(val) * weight; Mean::push(val, weight); }
	double sum2() const { return Sum2; }
	double variance() const { return (Sum2 - sqr(Sum) / Count) / (Count - 1); }
	double varianceMLE() const { return (Sum2 - sqr(Sum) / Count) / Count; }
	double stddeviation() const { return sqrt(variance()); }
	double stddeviationMLE() const { return sqrt(varianceMLE()); }
};

class Skewness : public Variance
{protected:
	double Sum3;
 public:
	Skewness() : Sum3(0) {}
	void push(double val, double weight = 1);
	double sum3() const { return Sum3; }
	double skewness() const;
	double skewnessBC() const { return skewness() * sqrt(Count * (Count - 1)) / (Count - 2); }
};

class Kurtosis : public Skewness
{protected:
	double Sum4;
 public:
	Kurtosis() : Sum4(0) {}
	void push(double val, double weight = 1) { Sum4 += sqr(sqr(val)) * weight; Skewness::push(val, weight); }
	double sum4() const { return Sum4; }
	double kurtosis() const;
	double kurtosisExcess() const { return kurtosis() - 3; }
	double kurtosisExcessBC() const;
	double kurtosisBC() const { return kurtosisExcessBC() + 3; }
};

#endif // MOMENT_H_
