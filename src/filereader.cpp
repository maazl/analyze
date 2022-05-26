#include "mathx.h"
#include <cmath>
#include "filereader.h"

using namespace std;


FileReader::FileReader(const char* filename)
:	FileName(filename)
,	In(filename, "r")
{}

void FileReader::SkipComments()
{	int n;
	do
	{	n = 0;
		ignore_result(fscanf(In, "#%*[^\n]\n%n", &n));
	} while (n);
}

bool FileReader::ReadLine(const unique_array<double>& dest)
{	if (!In)
		return false;
	SkipComments();
	if (feof(In))
	{	In = nullptr;
		return false;
	}
	for (auto& val : dest)
		if (fscanf(In, "%lg", &val) != 1)
			die(27, "Failed to read numeric data from file '%s' at byte %li: %s",
				FileName, ftell(In), ferror(In) ? strerror(errno) : feof(In) ? "unexpected end of file" : "no numeric input");
	ignore_result(fscanf(In, "%*[^\n]\n")); // discard remaining part of line.
	//fprintf(stderr, "ReadLine: %f\t%f\t%f\t%f\t%f\t%f\n", Input[0][0], Input[0][1], Input[0][2], Input[1][0], Input[1][1], Input[1][2]);
	return true;
}


void Interpolation::CalcInterpolation(double f)
{	for (unsigned i = 0; i < Result.size(); ++i)
		Result[i] = Input[0][i] * (1-f) + Input[1][i] * f;
	//fprintf(stderr, "Calc: %f\t%f\t%f\n", Result[0], Result[1], Result[2]);
}

Interpolation::Interpolation(unsigned count)
:	Input{ unique_array<double>(count + 1), unique_array<double>(count + 1) }
,	Result(count + 1)
{}

const unique_array<double>& Interpolation::Interpolate(double key)
{	double f = (key-Input[0][0]) / (Input[1][0]-Input[0][0]);
	//fprintf(stderr, "Interpolate: %f\t%f\t%f\t%f\t%i\n", Input[0][0], Input[1][0], key, f, Current);
	if (f <= 0.)
		return Input[0];
	if (f >= 1.)
		return Input[1];
	CalcInterpolation(f);
	return Result;
}

bool FileInterpolation::ReadLine()
{	if (FileReader::ReadLine(Feed()))
		return true;
	Feed();
	return false;
}

FileInterpolation::FileInterpolation(const char* filename, unsigned count)
:	Interpolation(count)
,	FileReader(filename)
{	if (!ReadLine() || !ReadLine())
		die(27, "Source file '%s' must contain at least 2 data lines.", FileName);
}

const unique_array<double>& FileInterpolation::Get(double key)
{	// Read enough data
	if (In)
		while (Last()[0] < key && ReadLine());
	return Interpolate(key);
}


void PolarFileInterpolation::CalcInterpolation(double f)
{	for (unsigned i = 1; i < Result.size(); i += 2)
	{	double abs0 = sqrt(sqr(Input[0][i]) + sqr(Input[0][i+1]));
		double arg0 = atan2(Input[0][i+1], Input[0][i]);
		double abs1 = sqrt(sqr(Input[1][i]) + sqr(Input[1][i+1]));
		double arg1 = atan2(Input[1][i+1], Input[1][i]);
		arg1 -= M_2PI * floor((arg1 - arg0)/M_2PI + .5);
		double abs = abs0 * (1-f) + abs1 * f;
		if (std::isnan(arg0)) arg0 = 0;
		if (std::isnan(arg1)) arg1 = 0;
		//double arg = (arg0 * (1-f) * abs0 + abs1 * f * abs1) / ((1-f) * abs0 + f * abs1);
		double arg = arg0 * (1-f) + arg1 * f;
		Result[i] = abs * cos(arg);
		Result[i+1] = abs * sin(arg);
	}
}
