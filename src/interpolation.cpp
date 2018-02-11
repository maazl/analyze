#include "interpolation.h"
#include "mathx.h"
#include <cmath>

using namespace std;


void Interpolation::CalcInterpolation(double f)
{	for (unsigned i = 1; i < Result.size(); ++i)
		Result[i] = Input[0][i] * (1-f) + Input[1][i] * f;
}

Interpolation::Interpolation(unsigned count)
:	Input{ unique_array<double>(count + 1), unique_array<double>(count + 1) }
,	Result(count + 1)
{}

const unique_array<double>& Interpolation::Interpolate(double key)
{	double f = (key-Input[0][0]) / (Input[0][0]-Input[1][0]);
	if (f <= 0.)
		return Input[0];
	if (f >= 1.)
		return Input[1];
	Result[0] = key;
	CalcInterpolation(f);
	return Result;
}

void FileInterpolation::SkipComments()
{	int n;
	do
	{	n = 0;
		ignore_result(fscanf(In, "#%*[^\n]\n%n", &n));
	} while (n);
}

void FileInterpolation::ReadLine()
{	if (!In)
		return;
	SkipComments();
	if (feof(In))
	{	In = nullptr;
		return;
	}
	auto& dest = Feed();
	for (auto& val : dest)
		if (fscanf(In, "%lg", &val) != 1)
			die(27, "Failed to read numeric data from file '%s' at byte %li: %s",
				FName, ftell(In), ferror(In) ? strerror(errno) : feof(In) ? "unexpected end of file" : "no numeric input");
	ignore_result(fscanf(In, "%*[^\n]\n")); // discard remaining part of line.
}

FileInterpolation::FileInterpolation(const char* filename, unsigned count)
:	Interpolation(count)
,	FName(filename)
,	In(filename, "r")
{	ReadLine();
	ReadLine();
	if (!In)
		die(27, "Source file '%s' must contain at least 2 data lines.", FName);
}

const unique_array<double>& FileInterpolation::Get(double key)
{	// Read enough data
	while (In && Last()[0] < key)
		ReadLine();
	return Interpolate(key);
}


void PolarFileInterpolation::CalcInterpolation(double f)
{	for (unsigned i = 1; i < Result.size(); i += 2)
	{	double abs0 = sqrt(sqr(Input[0][i]) + sqr(Input[0][i+1]));
		double arg0 = atan2(Input[0][i+1], Input[0][i]);
		double abs1 = sqrt(sqr(Input[1][i]) + sqr(Input[1][i+1]));
		double arg1 = atan2(Input[1][i+1], Input[1][i]);
		double abs = abs0 * (1-f) + abs1 * f;
		if (std::isnan(arg0)) arg0 = 0;
		if (std::isnan(arg1)) arg1 = 0;
		double arg = (arg0 * (1-f) * abs0 + abs1 * f * abs1) / (abs0 + abs1);
		Result[i] = abs * cos(arg);
		Result[i+1] = abs * sin(arg);
	}
}
