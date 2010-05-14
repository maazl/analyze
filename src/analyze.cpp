#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <386/builtin.h>
#include <stdarg.h>

#include <rfftw.h>
#include <complex>
#include <vector>
using namespace std;
typedef complex<double> Complex;
#include "pca.h"
#include "parser.h"

#define M_2PI (2.*M_PI)
#define M_3PI (3.*M_PI)
#ifndef M_PI_2
#define M_PI_2 (M_PI*.5)
#endif
#define M_PI_180 (M_PI/180.)
#define M_180_PI (180./M_PI)

#define N_MAX (65536*8)
#define CA_MAX 2
#define HA_MAX 5

//#define _cdecl // semms to make trouble

// avoid name clash with math.h
#define fmin fmin__
#define fmax fmax__


// data buffers
static short inbuffertmp[2*CA_MAX*N_MAX]; // Buffer for raw input
static float* inbuffer1 = NULL;  // Buffer for nominator input
static float* inbuffer2 = NULL;  // Buffer for denominator input
static float* ovrbuffer1 = NULL; // Buffer for overridden nominator
static float* ovrbuffer2 = NULL; // Buffer for overridden denominator
static float* outbuffer1 = NULL; // Buffer for FFT(inbuffer1)
static float* outbuffer2 = NULL; // Buffer for FFT(inbuffer2)
static float* ccbuffer1 = NULL;  // Buffer for cross correlation temporary data
static float* ccbuffer2 = NULL;  // Buffer for cross correlation of outbuffer
static float window[N_MAX];      // Buffer for window function
static int* harmonics = NULL;    // Buffer for harmonics dispatch table


static inline double sqr(double v)
{  return v*v;
}

static inline double todB(double f)
{  return 20*log10(f);
}

static inline double fromdB(double d)
{  return pow(10, d/20);
}

static double noiselvl_;

static const double minval = 1E-20;

static double getweight(double a1, double a2, double)
{  /*a1 -= noiselvl;
   a2 -= noiselvl;
   if (a1 <= 0)
      a1 = 0;
   if (a2 <= 0)
      a2 = 0;*/
   return noiselvl_/(1/sqr(a1+minval) + 1/sqr(a2+minval));
}

static double getweightD(double a1, double a2, double)
{  /*a1 -= noiselvl;
   a2 -= noiselvl;
   if (a1 <= 0)
      a1 = 0;
   if (a2 <= 0)
      a2 = 0;*/
   a1 += minval;
   double a1q = 1/sqr(a1);
   return noiselvl_/(a1q + sqrt(a1q + 1/sqr(a1 + a2))/(a2 + minval));
}

static double getconstweight(double, double, double)
{  return noiselvl_;
}

static double get1_fweight(double, double, double f)
{  return noiselvl_/f;
}

static void vectorscale(float* data, double factor, size_t len)
{  while (len--)
      *data++ *= factor;
}

static void vectoradd(float* dst, float* src, size_t len)
{  while (len--)
      *dst++ += *src++;
}

static void vektormul(float* dst, float* src, size_t len)
{  while (len--)
      *dst++ *= *src++;
}

// config
static float      gainadj[2]  = {1,1}; // gain {l, r}
static unsigned   N           = 8192;  // FFT length
static double     noiselvl    = 1;     // ?
static unsigned   winfn       = 0;     // window function: 0 = rectangle, 1 = Bartlett, 2 = Hanning, 3 = Hamming, 4 = Blackman, 5 = Blackman-Harris
static double     freq        = 48000; // sampling rate
static double     fmin        = -1;    // minimum freuqncy for FFT analysis
static double     fmax        = 1E99;  // minimum freuqncy for FFT analysis
static double     famin       = 1;     // ignore frquencies below famin for calculation of screen output
static double     famax       = 1E99;  // ignore frquencies above famax for calculation of screen output
static bool       writeraw    = false; // write raw data to file
static bool       writedata   = false; // write analysis data to file
static bool       writewindow = false; // write window function data to file
static unsigned   method      = 0;     // analysis method: 0 = none, 1 = PCA, 2 = FFT, 3 = PCA & FFT, 4 = XY
static unsigned   purgech     = 1;     // set the first FFT frequencies to 0
static unsigned   discsamp    = 0;     // skip the first samples
static bool       disctrail   = false; // comsume trailing samples after completion
static unsigned   addch       = 1;     // binsize in raw samples
static unsigned   addloop     = 1;     // add raw data before analysis # times
static bool       incremental = false; // incremental mode (add all raw data)
static double     rref        = 1;     // value of the reference resistor in impedance measurements
static unsigned   scalemode   = 1;     // l/r matrix decoder: 1 = L=l & R=r, 2 = L=r & R=l-r, 3 = L=r & R=l
static bool       stereo      = false; // Stereo aggregate mode (Toggle harmonics)
static unsigned   loops       = 1;     // number of analysis loops
static unsigned   lpause      = 10;    // number of loops between zero calibration parts
static unsigned   zeromode    = 0;     // zero calibration mode: 0 = none, 1 = read, 2 = generate, 3 = generatedelta, 4 = generate part 2, 5 = generatedelta part 2
static unsigned   gainmode    = 0;     // gain calibration mode: 0 = none, 1 = read, 2 = generate, 3 = generatedelta
static double     linphase    = 0;     // linear phase correction [s]
static bool       normalize   = false; // normalize L+R to 1. for impedance measurements
static unsigned   binsz       = 1;     // binsize in FFT channels
static double     fbinsc      = 0;     // logarithmic binsize: fmin/fmax = 1 + fbinsc
static double     f_inc       = 1;     // Absolute increment for harmonic table calculation
static double     f_log       = 0;     // Relative increment for harmonic table calculation
static unsigned   harmonic    = 0;     // analyze up to # harmonics
static bool       crosscorr   = false; // Calculate and remove time delay by cross correlation
static const char* datafile   = "data.dat"; // filename for analysis data
static const char* zerofile   = "zero.dat"; // file name for zero calibration data
static const char* zerodifffile="zeroD.dat";// file name for differential zero calibration data
static const char* gainfile   = "gain.dat"; // file name for gain calibration data
static const char* gaindifffile="gainD.dat";// file name for differential gain calibration data
static const char* rawfile    = "raw.dat";  // file name for raw data
static const char* windowfile = "window.dat";// file name for window data
static struct ovrwrt                   // overwrite channel with ...
{ const char* file;                    // ... file
  unsigned    column;                  // ... column in file
} overwrt[2] = {{NULL, 1}, {NULL, 1}};
static const char* infile     = "-";   // input file name
static const char* execcmd    = NULL;  // shell command to execute after analysis
static const char* plotcmd    = NULL;  // string to write to stdout after analysis
static double     (*weightfn)(double, double, double) = getweight;// weight function
// internal vars
static rfftw_plan P;    // FFT plan for forward transformation
static rfftw_plan PI;   // FFT plan for inverse transformation


static void createwindow(float* dst, int type, size_t len)
{  ++len;
   double sum = 0;
   float* win = dst;
   for (size_t i = 1; i < len; i++)
   {  double w;
      switch (type)
      {default:// rectangle
         w = 1;
         break;
       case 1: // Bartlett
         w = abs(i - len/2.);
         break;
       case 2: // Hanning
         w = .5 - .5*cos(2*M_PI*i/len);
         break;
       case 3: // Hamming
         w = .54 - .46*cos(2*M_PI*i/len);
         break;
       case 4: // Blackman
         w = .42 - .5*cos(2*M_PI*i/len) + .08*cos(4*M_PI*i/len);
         break;
       case 5: // Blackman-Harris
         w = .35875 - .48829*cos(2*M_PI*i/len) + .14128*cos(4*M_PI*i/len) - .01168*cos(6*M_PI*i/len);
      }
      sum += *win++ = w;
   }
   vectorscale(dst, len/sum, len);
}


static inline short fromraw(short v)
{  return _srotl(v, 8);
}

static int minmax[4];

static inline short storeminmax(short val, int* dst)
{  if (val < dst[0])
      dst[0] = val;
   if (val > dst[1])
      dst[1] = val;
   return val;
}

static void short2float(float* dst, const short* src, size_t len)
{  while (len)
   {  double d = 0;
      int i = addch;
      do d += storeminmax(fromraw(*src++), minmax);
       while (--i);
      *dst++ = d * gainadj[0];
      --len;
   }
}

static void short2float2(float* dst1, float* dst2, const short* src, size_t len)
{  while (len)
   {  double d1 = 0;
      double d2 = 0;
      int i = addch;
      do
      {  d1 += storeminmax(fromraw(src[0]), minmax);
         d2 += storeminmax(fromraw(src[1]), minmax+2);
         src += 2;
      } while (--i);
      *dst1++ = d1 * gainadj[0];
      *dst2++ = d2 * gainadj[1];
      --len;
   }
}

static void short2float2add(float* dst1, float* dst2, const short* src, size_t len)
{  while (len)
   {  double d1 = 0;
      double d2 = 0;
      int i = addch;
      do
      {  d1 += storeminmax(fromraw(src[0]), minmax);
         d2 += storeminmax(fromraw(src[1]), minmax+2);
         src += 2;
      } while (--i);
      *dst1++ += d1 * gainadj[0];
      *dst2++ += d2 * gainadj[1];
      --len;
   }
}

static void short2float2window(float* dst1, float* dst2, const short* src, size_t len)
{  const float* win = window;
   while (len)
   {  double d1 = 0;
      double d2 = 0;
      int i = addch;
      do
      {  d1 += storeminmax(fromraw(src[0]), minmax);
         d2 += storeminmax(fromraw(src[1]), minmax+2);
         src += 2;
      } while (--i);
      *dst1++ = d1 * *win * gainadj[0];
      *dst2++ = d2 * *win++ * gainadj[1];
      --len;
   }
}

static void short2floatD(float* dst1, float* dst2, const short* src, size_t len)
{  while (len)
   {  double d1 = 0;
      double d2 = 0;
      int i = addch;
      do
      {  d1 += storeminmax(fromraw(src[0]), minmax);
         d2 += storeminmax(fromraw(src[1]), minmax+2);
         src += 2;
      } while (--i);
      *dst2++ = d1 * gainadj[0]
       - (*dst1++ = d2 * gainadj[1]);
      --len;
   }
}

static void short2floatDadd(float* dst1, float* dst2, const short* src, size_t len)
{  while (len)
   {  double d1 = 0;
      double d2 = 0;
      int i = addch;
      do
      {  d1 += storeminmax(fromraw(src[0]), minmax);
         d2 += storeminmax(fromraw(src[1]), minmax+2);
         src += 2;
      } while (--i);
      d2 *= gainadj[1];
      *dst1++ += d2;
      *dst2++ += d1 * gainadj[0] - d2;
      --len;
   }
}

static void short2floatDwindow(float* dst1, float* dst2, const short* src, size_t len)
{  const float* win = window;
   while (len)
   {  double d1 = 0;
      double d2 = 0;
      int i = addch;
      do
      {  d1 += storeminmax(fromraw(src[0]), minmax);
         d2 += storeminmax(fromraw(src[1]), minmax+2);
         src += 2;
      } while (--i);
      *dst2++ = d1 * *win * gainadj[0]
       - (*dst1++ = d2 * *win * gainadj[1]);
      ++win;
      --len;
   }
}

static void applywindow(float* dst, size_t len)
{  const float* win = window;
   while (len)
      *dst++ *= *win++;
}

static void applywindowI(float* dst, size_t len)
{  const float* win = window;
   while (len)
      *dst++ /= *win++;
}

static inline double abs(double d1, double d2)
{  return sqrt(sqr(d1) + sqr(d2));
}

static void complex2polar(float* data, size_t len)
{  float *data2 = data + len;
   while (++data < --data2)
   {  register double phi = atan2(*data2, *data);
      *data = abs(*data, *data2);
      *data2 = phi;
   }
}

static void init()
{  minmax[0] = INT_MAX;
   minmax[1] = INT_MIN;
   minmax[2] = INT_MAX;
   minmax[3] = INT_MIN;
}

static void write1ch(FILE* out, const float* data, size_t len)
{  while (len--)
      fprintf(out, "%g\n", *data++);
}

static void write2ch(FILE* out, const short* data, size_t len)
{  while (len--)
   {  fprintf(out, "%i\t%i\n", fromraw(data[0]), fromraw(data[1]));
      data += 2;
   }
}

static void write2ch(FILE* out, const float* data, size_t len)
{  while (len--)
   {  fprintf(out, "%g\t%g\n", data[0], data[1]);
      data += 2;
   }
}

static void write2ch(FILE* out, const float* data1, const float* data2, size_t len)
{  while (len--)
      fprintf(out, "%g\t%g\n", *data1++, *data2++);
}

static void writepolar(FILE* out, const float* data, size_t len, double inc)
{  const float* data2 = data + len;
   // 1st line
   fprintf(out, "0\t%g\t0\n", *data++);
   len = 1;
   while (data < --data2)
      fprintf(out, "%g\t%g\t%g\n", len++*inc, *data++, *data2);
}

static void writecomplex(FILE* out, const Complex* data, size_t len)
{  while (len--)
   {  fprintf(out, "%14g\t%14g\t%14g\t%14g\n", data->real(), data->imag(), abs(*data), arg(*data)*M_180_PI);
      ++data;
}  }

static void write4complex(FILE* out, const Complex (* data)[4], size_t len)
{  while (len--)
   {  //Complex det = (*data)[0] * (*data)[3] - (*data)[1] * (*data)[2];
      fprintf(out, "%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\n"
       , (*data)[0].real(), (*data)[0].imag(), (*data)[1].real(), (*data)[1].imag(), (*data)[2].real(), (*data)[2].imag(), (*data)[3].real(), (*data)[3].imag()
       , abs((*data)[0]), arg((*data)[0])*M_180_PI, abs((*data)[1]), arg((*data)[1])*M_180_PI, abs((*data)[2]), arg((*data)[2])*M_180_PI, abs((*data)[3]), arg((*data)[3])*M_180_PI);
       //, det.real(), det.imag());
      ++data;
}  }

static void readcomplex(FILE* in, Complex* data, size_t len)
{  while (len--)
   {  double a,b;
      if (fscanf(in, "%lg%lg%*[^\n]", &a, &b) != 2)
         die(27, "Failed to read complex data (%i).", errno);
      //(stderr, "%g\t%g\n", a,b);
      *data++ = Complex(a,b);
}  }

static void read4complex(FILE* in, Complex (* data)[4], size_t len)
{  while (len--)
   {  double a,b,c,d,e,f,g,h;
      if (fscanf(in, "%lg%lg%lg%lg%lg%lg%lg%lg%*[^\n]", &a, &b, &c, &d, &e, &f, &g, &h) != 8)
         die(27, "Failed to read complex data (%i).", errno);
      //(stderr, "%g\t%g\n", a,b);
      (*data)[0] = Complex(a,b);
      (*data)[1] = Complex(c,d);
      (*data)[2] = Complex(e,f);
      (*data)[3] = Complex(g,h);
      ++data;
}  }

static void readfloat_2(FILE* in, unsigned column, size_t count, float* dest, size_t inc = 1)
{  while (count--)
   {  unsigned col = column;
      while (--col)
         fscanf(in, "%*s");
      if (fscanf(in, "%g%*[^\n]", dest) != 1)
         die(27, "Failed to raed column %u from data file.", column);
      dest += inc;
   }
}

/*static char* abbrev(const char* s, const char* token)
{  size_t l = strlen(s);
   return strnicmp(s, token, l) == 0 ? (char*)token + l : NULL;
}*/

static void readN(const char* s, unsigned* r)
{  bool ex;
   if (ex = *s == '^')
      ++s;
   readuint(s, r);
   if (ex)
      *r = 1 << *r;
}


static void dofft()
{  // forwardtransformation
   rfftw_one(P, inbuffer1, outbuffer1);
   rfftw_one(P, inbuffer2, outbuffer2);

   static const double minscale = 1E-15;
   if (purgech)
   {  // purge DC (meaningless without DC-coupling)
      outbuffer1[0] *= minscale;
      outbuffer2[0] *= minscale;
      for (int i = purgech; --i;)
      {  outbuffer1[i] *= minscale;
         outbuffer2[i] *= minscale;
         outbuffer1[N-i] *= minscale;
         outbuffer2[N-i] *= minscale;
   }  }
   // append constant zero to FFT result to simplify analysis
   outbuffer1[N] = 0;
   outbuffer2[N] = 0;
   
   // Calculate cross correlation to compensate for constant group delay.
   if (crosscorr)
   {  // calculate outbuffer1[] * conjugate(outbuffer2[])
      // f(0)
      ccbuffer1[0] = outbuffer1[0] * outbuffer2[0];
      // f(1..N/2-1)
      for (unsigned i = 1; i < N/2; ++i)
      { Complex c = Complex(outbuffer1[i], outbuffer1[N-i])
                  * Complex(outbuffer2[i], -outbuffer2[N-i]);
        ccbuffer1[i] = c.real();
        ccbuffer1[N-i] = c.imag();
      }
      // f(N/2)
      ccbuffer1[N/2] = outbuffer1[N/2] * outbuffer2[N/2];

      // do the cross correlation
      rfftw_one(PI, ccbuffer1, ccbuffer2);
      /*FILE* F = fopen("cc.dat", "w");
      write1ch(F, ccbuffer2, N);
      fclose(F);*/

      // use the result
      double phiinc = M_2PI / N;
      double asum = 0;
      double bsum = 0;
      double sum = 0;
      for (unsigned i = 0; i < N; ++i)
      { double amp = ccbuffer2[i];
        amp *= amp;
        asum += cos(phiinc * i) * amp;
        bsum += sin(phiinc * i) * amp;
        sum += amp;
      }
      asum /= sum;
      bsum /= sum;
      /*fprintf(stderr, "***** %12g %12g %12g %12g\n",
        sqrt(asum*asum+bsum*bsum), atan2(bsum, asum)*M_180_PI, asum, bsum);*/
  
      // calc linphase to compensate for the group delay
      linphase = atan2(bsum, asum) * N / freq;
   }
}

static double* wsums;
static Complex* gain;
static Complex* gainD;
static Complex (* zero)[4];
static Complex (* zeroD)[4];

// Calibration
static void docal(int bin, double f, Complex& U, Complex& I)
{  switch (gainmode)
   {case 1: // read
      U /= gain[bin];
      /*t = I * gainD[len];
      I -= U * gainD[len];
      U -= t;*/
      break;
    case 2: // write
    case 3:
      {  double weight = sqrt((*weightfn)(abs(U), abs(I), f));
         wsums[bin] += weight;
         U /= gain[bin];
         gainD[bin] += I/U * weight;
      }
   }
   if (normalize)
   {  Complex s = 1./(U + I);
      U *= s;
      I *= s;
   }
   if (zeromode & 1)
   {  Complex* cp = zero[bin];
      Complex t = U; // multiply (U,I) by (*cp)^(-1). The matrix inversion is easy because det(*cp) == 1.
      //Complex det = cp[0]*cp[3] - cp[1]*cp[2];
      U = (U * cp[3] - I * cp[1]);
      I = (-t * cp[2] + I * cp[0]);
      //U = (U * c22 - I * c12);
      //I = (-t * c21 + I * c11);
   }
   switch (zeromode)
   {case 2: // write part 1
    case 3:
      zeroD[bin][0] += U;
      zeroD[bin][2] += I;
      break;
    case 4: // write part 2
    case 5:
      zeroD[bin][1] += U;
      zeroD[bin][3] += I;
   }
}

// Phase unwrapper
// Adjusts phase by adding a multiple of 2pi so that the result
// is as close as possible to lph.
//   lph     last phase
//   phase   current calculated phase
//   returns unwrapped phase
static double unwrap(double lph, double phase)
{ return phase - M_2PI * floor((phase-lph) / M_2PI +.5);
}

class FFTbin
{public:
   enum    StoreRet
   {  BelowMin,  // frequency less than fmin
      AboveMax,  // frequency above fmax
      Ready,     // calculated values available
      Aggregated,// bin used for aggregation only
      Skip       // skip this bin because it is a harmonic
   };
 private:
   struct  aggentry
   {  double   f;
      double   Uabs; // Magnitude of nomiator
      double   Uarg; // Phase of nominator
      double   Iabs; // Magnitude of denominator
      double   Iarg; // Phase of denomiator 
      double   Zabs; // Magnitude of quotient
      double   Zarg; // Phase of quotient
      double   D;    // Group delay
      double   W;    // weight sum
      // internals
      double   lf;   // last frequency (for numerical derivative)
      double   lUarg;// last phase of nomiator
      double   lIarg;// last phase of denomiator
      double   lZarg;// last phase of quotient
      double   fnext;// next frequency for bin size
      unsigned binc; // number of bins accumulated
   };
 private:
   const double finc;

   aggentry  agg[2*HA_MAX+1];
   int       ch;
   aggentry* curagg;
   Complex   Zcache;

 public:
   FFTbin(double finc) : finc(finc)
   {  memset(agg, 0, sizeof agg);
   }
   StoreRet StoreBin(unsigned bin);

   void    PrintBin(FILE* dst) const;

   double   f() const { return curagg->f; } // frequency
   Complex  U() const { return polar(curagg->Uabs, curagg->Uarg); } // voltage, nominator or wanted signal
   double   Uabs() const { return curagg->Uabs; } // voltage magnitude
   double   Uarg() const { return curagg->Uarg; } // voltage phase
   Complex  I() const { return polar(curagg->Iabs, curagg->Iarg); } // current, denominator or reference signal
   double   Iabs() const { return curagg->Iabs; } // current magnitude
   double   Iarg() const { return curagg->Iarg; } // current phase
   Complex  Z() const { return Zcache; } // impedance, quotient or relative signal
   double   Zabs() const { return curagg->Zabs; } // impedance magnitude
   double   Zarg() const { return curagg->Zarg; } // impedance phase
   double   D() const { return curagg->D / M_2PI; } // group delay
   double   W() const { return curagg->W; } // weight
   int      h() const { return ch; }        // harmonic
 private:
};

FFTbin::StoreRet FFTbin::StoreBin(unsigned bin)
{  curagg = NULL;
   // frequency
   double f = bin * finc;
   ch = harmonics[bin];
   if ((unsigned)(abs(ch)-1) >= HA_MAX)
      return Skip;
   curagg = agg + ch + HA_MAX;
   unsigned base = bin;//ch != 0 ? bin/abs(ch) : bin; 
   // retrieve coefficients
   Complex U(outbuffer1[bin], bin && bin != N/2 ? outbuffer1[N-bin] : 0);
   Complex I(outbuffer2[base], bin && bin != N/2 ? outbuffer2[N-base] : 0);
   // calibration
   docal(bin, f, U, I);
   // phase correction
   U *= Complex(cos(linphase*f), sin(linphase*f));
   // calc Y
   Complex Z(U/I);
   // convert to polar
   double Uabs = abs(U);
   double Uarg = unwrap(curagg->lUarg, arg(U));
   double Iabs = abs(I);
   double Iarg = unwrap(curagg->lIarg, arg(I));
   double Zabs = abs(Z);
   double Zarg = unwrap(curagg->lZarg, arg(Z));
   // group delay
   double D = (Zarg - curagg->lZarg) / (f - curagg->lf);
   // store values for next point
   curagg->lUarg = Uarg;
   curagg->lIarg = Iarg;
   curagg->lZarg = Zarg; 
   curagg->lf    = f;

   // weight
   double w = (*weightfn)(Uabs, Iabs, f);
   if (curagg->binc == 0)
   {  // init
      curagg->f    = f * w;
      curagg->Uabs = Uabs * w;
      curagg->Uarg = Uarg * w;
      curagg->Iabs = Iabs * w;
      curagg->Iarg = Iarg * w;
      curagg->Zabs = Zabs * w;
      curagg->Zarg = Zarg * w;
      curagg->D    = D * w;
      curagg->W    = w;
      curagg->fnext = f * (1+fbinsc) - finc;
   } else
   {  curagg->f    += f * w;
      curagg->Uabs += Uabs * w;
      curagg->Uarg += Uarg * w;
      curagg->Iabs += Iabs * w;
      curagg->Iarg += Iarg * w;
      curagg->Zabs += Zabs * w;
      curagg->Zarg += Zarg * w;
      curagg->D    += D * w;
      curagg->W    += w;
   }
   ++curagg->binc;
   if (f < curagg->fnext)
      return Aggregated;
   w = curagg->W;
   curagg->f    /= w;
   curagg->Uabs /= w;
   curagg->Uarg /= w;
   curagg->Iabs /= w;
   curagg->Iarg /= w;
   curagg->Zabs /= w;
   curagg->Zarg /= w;
   curagg->D    /= w;
   curagg->binc = 0;
   Zcache = polar(curagg->Zabs, curagg->Zarg);
   /*if (curagg->f < fmin)
      return BelowMin;
   if (curagg->f > fmax)
      return AboveMax;*/
   return Ready;
}

void FFTbin::PrintBin(FILE* dst) const
{  fprintf(dst, "%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %6i\n",
   // f    |Hl|      phil               |Hr|      phir
      f(), Uabs(), Uarg()*M_180_PI, Iabs(), Iarg()*M_180_PI,
   // |Hl|/|Hr| phil-phir          re          im
      Zabs(), Zarg()*M_180_PI, Z().real(), Z().imag(),
   // weight delay harmonic
      W(), D(), h());
}


const ArgMap argmap[] = // must be sorted
{  {"ainc", (ArgFn)&setflag,    &incremental, true}
 , {"al",   (ArgFn)&readuint,   &addloop,     0}
 , {"bin",  (ArgFn)&readuint,   &binsz,       0}
 , {"bn" ,  (ArgFn)&readN,      &N,           0}
 , {"ca" ,  (ArgFn)&readuintdef,&addch,       2}
 , {"df" ,  (ArgFn)&readstring, &datafile,    0}
 , {"exec", (ArgFn)&readstring, &execcmd,     0}
 , {"famax",(ArgFn)&readdouble, &famax,       0}
 , {"famin",(ArgFn)&readdouble, &famin,       0}
 , {"fbin", (ArgFn)&readdouble, &fbinsc,      0}
 , {"finc", (ArgFn)&readdouble, &f_inc,       0}
 , {"flog", (ArgFn)&readdouble, &f_log,       0}
 , {"fmax", (ArgFn)&readdouble, &fmax,        0}
 , {"fmin", (ArgFn)&readdouble, &fmin,        0}
 , {"fq" ,  (ArgFn)&readdouble, &freq,        0}
 , {"g2f" , (ArgFn)&readstring, &gaindifffile,0}
 , {"gd" ,  (ArgFn)&setuint,    &gainmode,    3}
 , {"gf" ,  (ArgFn)&readstring, &gainfile,    0}
 , {"gg" ,  (ArgFn)&setuint,    &gainmode,    2}
 , {"gr" ,  (ArgFn)&setuint,    &gainmode,    1}
 , {"h/f",  (ArgFn)&setuint,    &weightfn,    (int)get1_fweight}
 , {"harm", (ArgFn)&readuint,   &harmonic,    0}
 , {"hd" ,  (ArgFn)&setuint,    &weightfn,    (int)getweightD}
 , {"he" ,  (ArgFn)&setuint,    &weightfn,    (int)getconstweight}
 , {"in" ,  (ArgFn)&readstring, &infile,      0}
 , {"ln" ,  (ArgFn)&readuint,   &loops,       1}
 , {"loop", (ArgFn)&setuint,    &loops,       INT_MAX}
 , {"lp" ,  (ArgFn)&readuint,   &lpause,      0}
 , {"lvl",  (ArgFn)&readdouble, &noiselvl,    0}
 , {"mfft", (ArgFn)&setbit,     &method,      1}
 , {"mpca", (ArgFn)&setbit,     &method,      2}
 , {"mst",  (ArgFn)&setflag,    &stereo,      true}
 , {"mxy",  (ArgFn)&setbit,     &method,      4}
 , {"olc",  (ArgFn)&readuint,   &overwrt[0].column, 0}
 , {"olf",  (ArgFn)&readstring, &overwrt[0].file,   0}
 , {"orc",  (ArgFn)&readuint,   &overwrt[1].column, 0}
 , {"orf",  (ArgFn)&readstring, &overwrt[1].file,   0}
 , {"pdc",  (ArgFn)&readuintdef,&purgech,     1}
 , {"phcc", (ArgFn)&setflag,    &crosscorr,   true}
 , {"phl",  (ArgFn)&readdouble, &linphase,    0}
 , {"plot", (ArgFn)&readstring, &plotcmd,     0}
 , {"psa",  (ArgFn)&readuintdef,&discsamp,    1}
 , {"pte",  (ArgFn)&readuintdef,&disctrail,   1}
 , {"rf" ,  (ArgFn)&readstring, &rawfile,     0}
 , {"rref", (ArgFn)&readdouble, &rref,        0}
 , {"scm",  (ArgFn)&readuint,   &scalemode,   0}
 , {"wd" ,  (ArgFn)&setflag,    &writedata,   true}
 , {"wf" ,  (ArgFn)&readstring, &windowfile,  0}
 , {"win",  (ArgFn)&readuintdef,&winfn,       2}
 , {"wr" ,  (ArgFn)&setflag,    &writeraw,    true}
 , {"ww" ,  (ArgFn)&setflag,    &writewindow, true}
 , {"z2f" , (ArgFn)&readstring, &zerodifffile,0}
 , {"zd" ,  (ArgFn)&setuint,    &zeromode,    3}
 , {"zf" ,  (ArgFn)&readstring, &zerofile,    0}
 , {"zg" ,  (ArgFn)&setuint,    &zeromode,    2}
 , {"zn" ,  (ArgFn)&setuint,    &normalize,   1}
 , {"zr" ,  (ArgFn)&setuint,    &zeromode,    1}
};
const size_t argmap_size = sizeof argmap / sizeof *argmap;

int main(int argc, char* argv[])
{  // parse cmdl
   while (--argc)
      parsearg(*++argv);

   if (N > N_MAX)
      die(32, "FFT Length too large.");
   if (N * addch > N_MAX * CA_MAX)
      die(32, "Input data length too large.");
   if (overwrt[0].file && overwrt[1].file)
      infile = NULL;

   // prepare some global vars
   noiselvl_ = 1/noiselvl;
   freq /= addch;
   linphase *= M_2PI;
   f_inc -= .5;
   f_log += 1;
   gainadj[0] /= 32767;
   gainadj[1] /= 32767;
   // allocate buffers
   if (method & 4)
   {  // reserve space for integrals and differentials too
      inbuffer1 = new float[3*N];
      inbuffer2 = new float[3*N];
      outbuffer1 = new float[3*N+1];
      outbuffer2 = new float[3*N+1];
   } else
   {  inbuffer1 = new float[N];
      inbuffer2 = new float[N];
      outbuffer1 = new float[N+1];
      outbuffer2 = new float[N+1];
   }
   if (crosscorr && (method & 1))
   {  ccbuffer1 = new float[N];
      ccbuffer2 = new float[N];
   }
   harmonics = new int[N/2+1];
   wsums = new double[N_MAX/2+1];
   gain = new Complex[N_MAX/2+1];
   gainD = new Complex[N_MAX/2+1];
   zero = new Complex[N_MAX/2+1][4];
   zeroD = new Complex[N_MAX/2+1][4];
   // create plan
   //float in[N], tout[N], power_spectrum[N/2+1];
   P = rfftw_create_plan(N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
   PI = rfftw_create_plan(N, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
   // create harmonics table
   {  memset(harmonics, 0, (N/2+1) * sizeof *harmonics);
      int sign = 1;
      for (unsigned i = (int)floor(fmin/freq*N +.5); i <= floor(fmax/freq*N +.5); ++i)
      {  if (i)
         {  for (unsigned j = 1; j <= harmonic && i*j <= N/2; ++j)
               if (harmonics[i*j])
                  goto next_f;
            for (unsigned j = 1; i*j <= N/2; ++j)
               harmonics[i*j] = j * sign;
         }
         if (stereo)
            sign = -sign;
         i = (int)floor(i * f_log + f_inc);
       next_f:;
      }
      /*FILE* f = fopen("harm.dat", "w");
      for (unsigned i = 0; i <= N/2; ++i)
         fprintf(f, "%12g %8i\n", i*freq/N, harmonics[i]);
      fclose(f);*/ 
   }

   createwindow(window, winfn, N);
   // write window data
   if (writewindow)
   {  FILE* tout;
      tout = fopen(windowfile, "wt");
      if (tout == NULL)
         die(21, "Failed to open %s for writing.", windowfile);
      write1ch(tout, window, N);
      fclose(tout);
   }

   FILE* in = NULL;
   if (infile)
   { if (strcmp(infile, "-") == 0)
     {  // streaming
        in = stdin;
        _fsetmode(in, "b");
     } else
     {  in = fopen(infile, "rb");
        if (in == NULL)
           die(20, "Failed to open input file.");
     }
     // discard first samples
     fread(inbuffertmp, sizeof *inbuffertmp, discsamp, in);
   }
   
   if (overwrt[0].file)
   {  FILE* fin = fopen(overwrt[0].file, "r");
      if (fin == NULL)
         die(20, "Failed to open %s for reading.", overwrt[0].file);
      ovrbuffer1 = new float[N];
      readfloat_2(fin, overwrt[0].column, N, ovrbuffer1);
      fclose(fin);
   }
   if (overwrt[1].file)
   {  FILE* fin = fopen(overwrt[1].file, "r");
      if (fin == NULL)
         die(20, "Failed to open %s for reading.", overwrt[1].file);
      ovrbuffer2 = new float[N];
      readfloat_2(fin, overwrt[0].column, N, ovrbuffer2);
      fclose(fin);
   }
   
   memset(inbuffer1, 0, sizeof inbuffer1); // init with 0 because of incremental mode
   memset(inbuffer2, 0, sizeof inbuffer2);
   memset(wsums, 0, sizeof wsums);
   memset(wsums, 0, sizeof wsums);
   memset(gainD, 0, sizeof gain);
   memset(zeroD, 0, sizeof zeroD);
   // prepare gainmode
   FILE* fz;
   switch (gainmode)
   {case 1: // read
    case 3:
      fz = fopen(gainfile, "r");
      if (fz == NULL)
         die(20, "Failed to open %s for reading.", gainfile);
      readcomplex(fz, gain, N/2+1);
      fclose(fz);
   }
 restart_zero:
   // prepare zeromode
   switch (zeromode)
   {case 1: // read
    case 3:
      FILE* fz;
      fz = fopen(zerofile, "r");
      if (fz == NULL)
         die(20, "Failed to open %s for reading.", zerofile);
      read4complex(fz, zero, N/2);
      fclose(fz);
   }

   // operation loop
   unsigned loop = loops;
   unsigned addloops = 0;
   do
   {  
      if (in)
      {  if (fread(inbuffertmp, 2*sizeof *inbuffertmp * addch, N, in) != N)
            die(27, "Failed to read data from input.");

         // write raw data
         if (writeraw)
         {  FILE* tout = fopen(rawfile, "wt");
            if (tout == NULL)
               die(21, "Failed to open %s for writing.", rawfile);
            write2ch(tout, inbuffertmp, N * addch);
            fclose(tout);
         }

         // reset min/max
         init();

         switch (scalemode)
         {case 1:
            if (incremental || addloops)
               short2float2add(inbuffer1, inbuffer2, inbuffertmp, N);
             else if (winfn && addloop == 1)
               short2float2window(inbuffer1, inbuffer2, inbuffertmp, N);
             else
               short2float2(inbuffer1, inbuffer2, inbuffertmp, N);
            break;
          case 2:
            if (incremental || addloops)
               short2floatDadd(inbuffer1, inbuffer2, inbuffertmp, N);
             else if (winfn && addloop == 1)
               short2floatDwindow(inbuffer1, inbuffer2, inbuffertmp, N);
             else
               short2floatD(inbuffer1, inbuffer2, inbuffertmp, N);
            break;
          case 3:
            if (incremental || addloops)
               short2float2add(inbuffer2, inbuffer1, inbuffertmp, N);
             else if (winfn && addloop == 1)
               short2float2window(inbuffer2, inbuffer1, inbuffertmp, N);
             else
               short2float2(inbuffer2, inbuffer1, inbuffertmp, N);
            break;
          default:
            die(32, "Invalid scalmode");
         }
         // write raw status
         fprintf(stderr, "\nmin:\t%i\t%i\nmax:\t%i\t%i\n", minmax[0], minmax[2], minmax[1], minmax[3]);
      }
      
      if (incremental || addloops)
      {  if (ovrbuffer1)
            vectoradd(inbuffer1, ovrbuffer1, N);
         if (ovrbuffer2)
            vectoradd(inbuffer2, ovrbuffer2, N);
      } else
      {  if (ovrbuffer1)
            memcpy(inbuffer1, ovrbuffer1, N * sizeof *inbuffer1);
         if (ovrbuffer2)
            memcpy(inbuffer2, ovrbuffer2, N * sizeof *inbuffer2);
      }

      /*{  FILE* f = fopen("ov.dat", "w");
         for (size_t i = 0; i < N; ++i)
            fprintf(f, "%12g %12g\n", inbuffer1[i], inbuffer2[i]);
         fclose(f);
      }*/
      
      if (++addloops < addloop)
         continue; // add more data
      if (winfn && (incremental || addloops))
      {  // optimization: apply window function later
         applywindow(inbuffer1, N);
         applywindow(inbuffer2, N);
      }
      addloops = 0;

      // write raw data
      /*if (writeraw)
      {  tout = fopen(rawfile, "wt");
         write2ch(tout, inbuffer1, inbuffer2, N);
         fclose(tout);
      }*/

      switch (method)
      {default:
         die(34, "Invalid combination of measurement modes (%i), e.g. FFT and XY.", method);

       case 2: // PCA analysis
         {  PCA<5> pca;
            double data[6];
            float* U = inbuffer1 +1;
            float* I = inbuffer2 +1;
            data[2] = 1;   // offset
            data[3] = 0;   // integral
            data[4] = 0;   // linear
            data[5] = 0;   // differential
            unsigned i = N -3;
            do
            {  data[0] = U[0] + U[1];
               data[1] = I[0] + I[1];
               data[3] += I[-1] + I[0];
               data[5] = I[-1] + I[0] - I[1]  - I[2];
               pca.Store(*(double (*)[5])&data, (*weightfn)(data[0], data[1], i));
               data[4]++;
               U += 2;
               I += 2;
               i -= 2;
            } while (i > 0);

            Vektor<4> res = pca.Result()*rref;
            printf("\nPCA: %12g %12g %12g %12g %12g %12g\n", res[0], res[1], 2./freq/res[2], res[3], 1./2*freq*res[2]/M_2PI/res[0], freq*res[4]/2);
         }
         break;

       case 1: // FFT
         {  dofft();

            vectorscale(outbuffer1, sqrt(1./N)/addch, N);
            vectorscale(outbuffer2, sqrt(1./N)/addch, N);

            // calc some sums
            double wsum = 0;
            int nsum = 0;
            double Rsum = 0;
            double R2sum = 0;
            double L2sum = 0;
            //double LCsum = 0; ==> -2 * wsum
            double C2sum = 0;
            double Lsum = 0;
            double Csum = 0;
            double d2sum = 0;
            //double RWsum = 0;
            // write data
            FILE* tout = NULL;
            if (writedata)
            {  tout = fopen(datafile, "wt");
               if (tout == NULL)
                  die(21, "Failed to create %s.", datafile);
            }

            FFTbin calc(freq/N);

            for (size_t len = 0; len <= N/2; ++len)
            {  // do calculations and aggregations
               switch (calc.StoreBin(len))
               {case FFTbin::AboveMax:
                  // write
                  if (tout && abs(calc.h()) > 1)
                     calc.PrintBin(tout);
                default:
                //case FFTbin::BelowMin:
                //case FFTbin::Aggregated:
                //case FFTbin::Skip:
                  continue;
                case FFTbin::Ready:
                  // write
                  if (tout)
                     calc.PrintBin(tout);
               }

               if (calc.f()/calc.h() < famin && calc.f()/calc.h() >= famax)
                  continue;
               // average
               ++nsum;
               wsum += calc.W();
               // resistivity
               Rsum += calc.W() * calc.Z().real();
               R2sum += calc.W() * sqr(calc.Z().real());
               // L & C
               L2sum += calc.W() * sqr(calc.f());
               // LCsum += weight; == wsum
               C2sum += calc.W() / sqr(calc.f());
               Lsum += calc.W() * calc.Z().imag() * calc.f();
               Csum += calc.W() * calc.Z().imag() / calc.f();
               d2sum = calc.W() * sqr(calc.Z().imag());
            }
            if (tout)
               fclose(tout);

            // calculate summary
            double R = Rsum / wsum;
            double RE = sqrt((R2sum/wsum - sqr(R))/(nsum -1));
            double sLC = 1/(sqr(wsum) - C2sum*L2sum);
            double L = (C2sum*Lsum - Csum*wsum) * sLC;
            double C = (Csum*L2sum - wsum*Lsum) * sLC;
            double LE = sqrt((C2sum*d2sum - sqr(Csum) - (d2sum*sqr(wsum) - 2*Csum*wsum*Lsum + C2sum*sqr(Lsum))/L2sum) * sLC / (nsum -2));
            double CE = sqrt((L2sum*d2sum - sqr(Lsum) - (d2sum*sqr(wsum) - 2*Csum*wsum*Lsum + L2sum*sqr(Csum))/C2sum) * sLC / (nsum -2));

            // write summary
            fprintf(stderr, "\n%6i %12g %12g %12g %12g %12g %12g %12g %12g\n", nsum, wsum, Rsum, R2sum, Csum, C2sum, Lsum, L2sum, sLC);
            fprintf(stderr, "\nreal (R)     \t%12g Ò %g\n"
                            "imaginary (C)\t%12g Ò %g\n"
                            "imaginary (L)\t%12g Ò %g\n"
             , rref * R, rref * RE
             , 1/(rref*C*M_2PI), 1 / (CE * rref*M_2PI)
             , rref*L/M_2PI, LE * rref / M_2PI );
            if (crosscorr)
              fprintf(stderr, "delay        \t%12g\n", linphase/M_2PI);

         }
         break;

       case 3: // FFT and then PCA
         {  // FFT
            dofft();

            vectorscale(outbuffer1, sqrt(1./N)/addch, N);
            vectorscale(outbuffer2, sqrt(1./N)/addch, N);

            // write data
            FILE* tout = NULL;
            if (writedata)
            {  tout = fopen(datafile, "wt");
               if (tout == NULL)
                  die(21, "Failed to create %s.", datafile);
            }

            PCA<2> pcaRe;
            PCA<3> pcaIm;
            double PCAdataRe[2];
            double PCAdataIm[3];
            // some values are const
            PCAdataRe[1] = 1;
            //PCAdataIm[1] = 1;

            FFTbin calc(freq/N);

            // 1st line
            for (size_t len = 0; len <= N/2; ++len)
            {  // do calculations and aggregations
               switch (calc.StoreBin(len))
               {case FFTbin::AboveMax:
                  // write
                  if (tout && calc.h() > 1)
                     calc.PrintBin(tout);
                default:
                //case FFTbin::BelowMin:
                //case FFTbin::Aggregated:
                //case FFTbin::Skip:
                  continue;
                case FFTbin::Ready:
                  // write
                  if (tout)
                     calc.PrintBin(tout);
               }

               if (calc.f()/calc.h() < famin && calc.f()/calc.h() >= famax)
                  continue;
               // component analysis
               PCAdataRe[0] = calc.Z().real();
               //PCAdataRe[2] = 1/af;
               //PCAdataRe[3] = f;
               PCAdataIm[0] = calc.Z().imag(); // fit imaginary part in conductivity
               PCAdataIm[1] = 1/calc.f();
               PCAdataIm[2] = calc.f();
               //PCAdataIm[3] = 1/af;
               //printf("Re: %12g %12g %12g %12g\n", PCAdataRe[0], PCAdataRe[1], PCAdataRe[2], weight);

               // add values
               pcaRe.Store(PCAdataRe, calc.W());
               pcaIm.Store(PCAdataIm, calc.W());
            }
            if (tout)
               fclose(tout);

            // calculate summary
            Vektor<1> resRe = pcaRe.Result();
            Vektor<2> resIm = pcaIm.Result();

            fprintf(stderr, "resRe: %12g %12g %12g\n", resRe[0], resRe[1], resRe[2]);
            //printf("resRe: %12g %12g %12g %12g\n", resRe[0], resRe[1], resRe[2], resRe[3]);
            fprintf(stderr, "resIm: %12g %12g %12g %12g\n", resIm[0], resIm[1], resIm[2], resIm[3]);

            double R0 = rref*resRe[0];
            //double R1_f = rref*resRe[1];
            double C0 = -1 / M_2PI / resIm[0] / rref;
            double L0 = resIm[1] / M_2PI * rref;
            //double C1f = -resIm[0] / M_2PI / rref;

            // write summary
            fprintf(stderr, "\nreal: R [Ohm]     \t%12g\n"
//                              "      R/f [Ohm/Hz]\t%12g\t%12g @100Hz\n"
                              "imaginary: C [ÊF] \t%12g\n"
//                              "      Cf [F Hz]   \t%12g\t%12g @100Hz\n"
                              "imaginary: L [ÊH] \t%12g\n"
//             , R0, R1_f, R1_f/100, C0, C1f, C1f*100);
               , R0, C0*1E6, L0*1E6);
            if (crosscorr)
              fprintf(stderr, "delay        \t%12g\n", linphase/M_2PI);
         }
         break;

       case 4: // XY mode
         {  // we need to do an FFT here, at least for the calibration
            dofft();
            // normalize (including retransformation)
            vectorscale(outbuffer1, sqrt(1./N/N)/addch, N);
            vectorscale(outbuffer2, sqrt(1./N/N)/addch, N);

            const double inc = freq/N;
            // U(f)
            float* a1 = outbuffer1;
            float* a2 = outbuffer2;
            float* b1 = a1 + N;
            float* b2 = a2 + N;
            *b1 = 0; // well, somewhat easier this way
            *b2 = 0;
            for (int len = 0; a1 < b1; len += harmonic, a1 += harmonic, a2 += harmonic, b1 -= harmonic, b2 -= harmonic)
            {  // calc
               double f = len*inc;
               Complex U(*a1, *b1);
               Complex I(*a2, *b2);
               // calibration
               docal(len, f, U, I);
               // store data
               *a1 = U.real();
               *b1 = U.imag();
               *a2 = I.real();
               *b2 = I.imag();
               // The integrals and differentials are calculated in the frequency domain.
               // This is at the expence of approximate a factor 2 computing time, since we have to do
               // one forward and one backward transformation for the zero compensation anyway.
               // The advantage is that we do not have to deal with the 1/2 time slot linear phse shift
               // of the numeric integration/differentiation in the time domain.
               Complex di(0, f); // differential operator
               // store integrals
               Complex C = U/di;
               a1[N] = C.real();
               b1[N] = C.imag();
               C = I/di;
               a2[N] = C.real();
               b2[N] = C.imag();
               // store differentials
               C = U*di;
               a1[2*N] = C.real();
               b1[2*N] = C.imag();
               C = I*di;
               a2[2*N] = C.real();
               b2[2*N] = C.imag();
            }
            // clear DC and nyquist frequencies of integral and diferential, since they do no allow 90¯ phase shift.
            outbuffer1[N] = outbuffer2[N] = 0;
            outbuffer1[2*N] = outbuffer2[2*N] = 0;
            outbuffer1[N+N/2] = outbuffer2[N+N/2] = 0;
            outbuffer1[2*N+N/2] = outbuffer2[2*N+N/2] = 0;

            // now do the inverse transform to get the corrected data back.
            rfftw_one(PI, outbuffer1, inbuffer1);
            rfftw_one(PI, outbuffer1 + N, inbuffer1 + N);
            rfftw_one(PI, outbuffer1 + 2*N, inbuffer1 + 2*N);
            rfftw_one(PI, outbuffer2, inbuffer2);
            rfftw_one(PI, outbuffer2 + N, inbuffer2 + N);
            rfftw_one(PI, outbuffer2 + 2*N, inbuffer2 + 2*N);

            // write data
            if (writedata)
            {  FILE* tout = fopen(datafile, "wt");
               if (tout == NULL)
                  die(21, "Failed to create %s.", datafile);

               const float* Up = inbuffer1;
               const float* Ip = inbuffer2;
               for (unsigned len = 0; len < N; ++len, ++Up, ++Ip)
                  fprintf(tout, "%12g %12g %12g %12g %12g %12g %12g\n",
                  // t         U    I    INT U  INT I  D U      D I
                     len/freq*harmonic, *Up, *Ip, Up[N], Ip[N], Up[2*N], Ip[2*N]);
               fclose(tout);
            }
         }
      }

      if (execcmd)
         system(execcmd);
      if (plotcmd)
      {  // for gnuplot!
         puts(plotcmd);
         fflush(stdout);
      }

      // undo window function because of incremental mode
      // TODO: this causes a loss of precision!
      if (winfn && incremental)
      {  applywindowI(inbuffer1, N);
         applywindowI(inbuffer2, N);
      }

   } while (--loop);

   if (disctrail)
      fread(inbuffertmp, 2*sizeof *inbuffertmp * addch, N, in);

   switch (gainmode)
   {case 2: // write
    case 3:
      double* wp = wsums + N/2+1;
      for (Complex* cp = gainD + N/2+1; --cp >= gainD; )
         *cp /= *--wp;  // scale 2 average
      FILE* fz;
      fz = fopen(gainmode == 3 ? gaindifffile : gainfile, "wb");
      if (fz == NULL)
         die(21, "Failed to open %s.", gainmode == 3 ? gaindifffile : gainfile);
      writecomplex(fz, gainD, N/2+1);
      fclose(fz);
   }
   switch (zeromode)
   {case 2: // write
    case 3:
      zeromode += 2;
      puts("Zeromode calibration part one is now complete.\n"
           "Part 2 will start at the end of the contdown.\7");
      for (int loop = lpause; loop;)
      {  printf("\r%u ", loop);
         if (fread(inbuffertmp, 2*sizeof *inbuffertmp * addch, N, in) != N)
            die(27, "failed to read data");
         --loop;
      }
      puts("\nNow at part 2.");
      goto restart_zero;
    case 4: // part 2
    case 5:
      for (Complex (* cp)[4] = zeroD + N/2+1; --cp >= zeroD; )
      {  // scale to fit det *cp == 1
         Complex det = 1. / sqrt((*cp)[0] * (*cp)[3] - (*cp)[1] * (*cp)[2]);
         if ((((*cp)[0] + (*cp)[3]) * det).real() < 0)
            det = -det;
         (*cp)[0] *= det;
         (*cp)[1] *= det;
         (*cp)[2] *= det;
         (*cp)[3] *= det;
      }
      FILE* fz;
      fz = fopen(zeromode == 5 ? zerodifffile : zerofile, "wb");
      if (fz == NULL)
         die(21, "Failed to open %s.", zeromode == 5 ? zerodifffile : zerofile);
      write4complex(fz, zeroD, N/2);
      fclose(fz);
   }

   puts("completed.");
   // read until EOF
   if (disctrail)
   {  //do
      //   puts("trail");
      while (fread(inbuffertmp, 2*sizeof *inbuffertmp * addch, N, in) > 0);
      //puts("fin");
   }

   // close stdin to signal playrec
   //fclose(in);

   return 0;
}

