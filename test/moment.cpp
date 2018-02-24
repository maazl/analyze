#include "../src/moment.h"
#include <assert.h>


int main()
{
	{	Kurtosis s;
		s.push(1, 2);
		s.push(-1, 2);
		assert(s.mean() == 0);
		assert(fabs(s.stddeviation() - 1.1547) < .0001);
		assert(s.skewness() == 0);
		assert(fabs(s.kurtosisExcess() + 2) < .0001);
	}

	{	Kurtosis s;
		s.push(2, 2);
		s.push(0, 2);
		assert(s.mean() == 1);
		assert(fabs(s.stddeviation() - 1.1547) < .0001);
		assert(s.skewness() == 0);
		assert(fabs(s.kurtosisExcess() + 2) < .0001);
	}

	{	Kurtosis s;
		s.push(-1, 2);
		s.push(2, 1);
		assert(s.mean() == 0);
		assert(fabs(s.stddeviation() - 1.73205) < .0001);
		assert(fabs(s.skewness() - M_SQRT1_2) < .0001);
		assert(fabs(s.kurtosisExcess() + 1.5) < .0001);
	}

	return 0;
}
