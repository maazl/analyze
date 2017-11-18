#ifndef MATHX_H_
#define MATHX_H_

#include <math.h>
#include <stdint.h>
#include <stdlib.h>


inline double sqr(double v)
{	return v * v;
}
inline unsigned sqr(int v)
{	return v * v;
}
inline uint64_t sqr(int64_t v)
{	return v * v;
}
inline double cb(double v)
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


#endif // MATHX_H_
