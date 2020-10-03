#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "unique_array.h"
#include "utils.h"
#include <cstdio>
#include <array>


/// Class to read numeric data from a text file.
/// Lines starting with # are ignored.
class FileReader
{public:
	/// File name @remarks mainly for diagnostic messages.
	const char* const FileName;
 protected:
	/// Input file handle or \c nullptr if EOF.
	FILEguard In;
 protected:
	/// Skip comment lines.
	void SkipComments();
 public:
	FileReader(const char* filename);
	/// Read a line of data from the input stream \see In.
	/// @param dest read data into this array.
	/// The number of elements to read is taken from the array size.
	/// Additional elements are ignored.
	/// @return false = EOF
	/// @post If an EOF condition occurs \c In is reset.
	bool ReadLine(const unique_array<double>& dest);
};

/// Helper class to interpolate numeric values from rows of data using the first column as index.
class Interpolation
{protected:
	/// The last two rows read from file.
	/// @remarks This is a (small) ring buffer.
	const std::array<const unique_array<double>, 2> Input;
	/// Result row returned by \see Interpolate.
	const unique_array<double> Result;
 private:
	/// Current index into \see Input with the \e last row read.
	bool Current = false;
 protected:
	/// @brief Do the interpolation from \see Input into \see Result.
	/// @details The default implementation does a linear interpolation and \e no extrapolation.
	/// It just returns the lower or the higher value in doubt.
	/// @param f Interpolation factor, normally in the range (0..1).
	/// Closer to 0.0 means closer to \c Input[0], closer to 1.0 means closer to \c Input[1].
	/// Values < 0 or > 1 request \e extrapolation.
	/// @pre The function must not be called before at least 2 rows have been filled with \see Feed().
	virtual void CalcInterpolation(double f);
 public:
	/// Create an interpolation.
	/// @param count Number of columns to interpolate \e excluding the first key column.
	Interpolation(unsigned count);
	virtual ~Interpolation() {}
	/// @brief Add a new row of data.
	/// @return Array with the data row to be filled.
	/// @details When the function returns you have to fill \e all data values in the returned array.
	/// No other member function should be called unless this is completed.
	/// The first value is the key column and must be strictly monotonic on consecutive calls.
	const unique_array<double>& Feed() { return Input[Current ^= true]; }
	/// Return the last row read. This is equivalent to the last value returned by \see Feed().
	/// @pre The function must not be called before at least one row have been filled with \see Feed().
	const unique_array<double>& Last() const { return Input[Current]; }
	/// Return the last row read. This is equivalent to the one before last value returned by \see Feed().
	/// @pre The function must not be called before at least two rows have been filled with \see Feed().
	const unique_array<double>& Prev() const { return Input[!Current]; }
	/// Get interpolated column values.
	/// @param key Key value to be used
	virtual const unique_array<double>& Interpolate(double key);
};

/// @brief Helper class to interpolate numeric values in a file using the first column as index.
/// @details The class does not read the entire file into memory.
/// Instead it uses stream processing to keep only needed data.
class FileInterpolation : protected Interpolation, protected FileReader
{protected:
	virtual bool ReadLine();
 public:
	/// Create an interpolation from filename.
	/// @param filename File name.
	/// @param count Number of columns to interpolate \e excluding the first key column.
	/// The file might contain additional columns. They are ignored.
	FileInterpolation(const char* filename, unsigned count);
	/// Get interpolated values.
	/// @param key Key value used for interpolation.
	/// The passed value must be strictly monotonic,
	/// i.e. larger than on any previous call to this function for this instance.
	/// @details The function automatically reads more data from the input stream as needed.
	const unique_array<double>& Get(double key);
};

/// Variant of \see FileInterpolation for complex numbers that does the interpolation in polar coordinates.
class PolarFileInterpolation : public FileInterpolation
{protected:
	virtual void CalcInterpolation(double f);
 public:
	/// Create an interpolation from filename.
	/// @param filename File name.
	/// @param count Number of complex number columns to interpolate \e excluding the first key column.
	/// Each complex number must span over two columns with the real component first and the imaginary part second.
	/// The key column is always real. So in fact \code 2 * count + 1 \endcode columns are read.
	/// The file might contain additional columns. They are ignored.
	PolarFileInterpolation(const char* filename, unsigned count) : FileInterpolation(filename, count << 1) {}
};


#endif // INTERPOLATION_H_
