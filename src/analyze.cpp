#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <386/builtin.h>
#include <stdarg.h>

#include <rfftw.h>
#include <complex>
using namespace std;
typedef complex<double> Complex;
#include "pca.h"

#define M_2PI (2.*M_PI)
#define M_3PI (3.*M_PI)
#ifndef M_PI_2
#define M_PI_2 (M_PI*.5)
#endif
#define M_PI_180 (M_PI/180.)
#define M_180_PI (180./M_PI)

#define N_MAX (65536*8)
#define CA_MAX 2

//#define _cdecl // semms to make trouble

// avoid name clash with math.h
#define fmin fmin__
#define fmax fmax__


// data buffers
static short inbuffertmp[2*CA_MAX*N_MAX];
static float* inbuffer1 = NULL;
static float* inbuffer2 = NULL;
static float* outbuffer1 = NULL;
static float* outbuffer2 = NULL;
static float window[N_MAX];


void die(const char* msg, ...)
{  va_list va;
   va_start(va, msg);
   vfprintf(stderr, msg, va);
   va_end(va);
   fputc('\n', stderr);
   exit(1);
}


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
static int        N           = 8192;  // FFT length
static double     noiselvl    = 1;     // ?
static int        winfn       = 0;     // window function: 0 = rectangle, 1 = Bartlett, 2 = Hanning, 3 = Hamming, 4 = Blackman, 5 = Blackman-Harris
static double     freq        = 44100; // sampling rate
static double     fmin        = -1;    // minimum freuqncy for FFT analysis
static double     fmax        = 1E99;  // minimum freuqncy for FFT analysis
static double     famin       = 1;     // ignore frquencies below famin for calculation of screen output
static double     famax       = 1E99;  // ignore frquencies above famax for calculation of screen output
static bool       writeraw    = false; // write raw data to file
static bool       writedata   = false; // write analysis data to file
static bool       writewindow = false; // write window function data to file
static int        method      = 0;     // analysis method: 0 = none, 1 = PCA, 2 = FFT, 3 = PCA & FFT, 4 = XY
static int        purgech     = 0;     // set the first FFT frequencies to 0
static int        discsamp    = 0;     // skip the first samples
static bool       disctrail   = false; // comsume trailing samples after completion
static int        addch       = 1;     // binsize in raw samples
static int        addloop     = 1;     // add raw data before analysis # times
static bool       incremental = false; // incremental mode (add all raw data)
static double     rref        = 1;     // value of the reference resistor in impedance measurements
static int        scalemode   = 1;     // l/r matrix decoder: 1 = L=l & R=r, 2 = L=r & R=l-r, 3 = L=r & R=l
static int        loops       = 1;     // number of analysis loops
static int        lpause      = 10;    // number of loops between zero calibration parts
static int        zeromode    = 0;     // zero calibration mode: 0 = none, 1 = read, 2 = generate, 3 = generatedelta, 4 = generate part 2, 5 = generatedelta part 2
static int        gainmode    = 0;     // gain calibration mode: 0 = none, 1 = read, 2 = generate, 3 = generatedelta
static double     linphase    = 0;     // lienar phase correction [s]
static bool       nophase     = false; // purge any phase information
static bool       normalize   = false; // normalize L+R to 1. for impedance measurements
static int        binsz       = 1;     // binsize in FFT channels
static double     fbinsc      = 0;     // logarithmic binsize: fmin/fmax = 1 + fbinsc
static int        harmonic    = 1;     // analyze only harmonics of this base frequency (FFT)
static const char* datafile   = "data.dat";  // filename for analysis data
static const char* zerofile   = "zero.dat";  // file name for zero calibration data
static const char* zerodifffile = "zeroD.dat";// file name for differential zero calibration data
static const char* gainfile   = "gain.dat";  // file name for gain calibration data
static const char* gaindifffile = "gainD.dat";// file name for differential gain calibration data
static const char* rawfile    = "raw.dat";   // file name for raw data
static const char* windowfile = "window.dat";// file name for window data
static const char* infile     = NULL;  // input file name
static const char* execcmd    = NULL;  // shell command to execute after analysis
static const char* plotcmd    = NULL;  // string to write to stdout after analysis
static double     (*weightfn)(double, double, double) = getweight;// weight function
// internal vars
static int        Nh;   // FFT length after harmonic extraction
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
      do d += storeminmax(fromraw(*src++), minmax) * gainadj[0];
       while (--i);
      *dst++ = d;
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

void init()
{  minmax[0] = INT_MAX;
   minmax[1] = INT_MIN;
   minmax[2] = INT_MAX;
   minmax[3] = INT_MIN;
}

void write1ch(FILE* out, const float* data, size_t len)
{  while (len--)
      fprintf(out, "%g\n", *data++);
}

void write2ch(FILE* out, const short* data, size_t len)
{  while (len--)
   {  fprintf(out, "%i\t%i\n", fromraw(data[0]), fromraw(data[1]));
      data += 2;
   }
}

void write2ch(FILE* out, const float* data, size_t len)
{  while (len--)
   {  fprintf(out, "%g\t%g\n", data[0], data[1]);
      data += 2;
   }
}

void write2ch(FILE* out, const float* data1, const float* data2, size_t len)
{  while (len--)
      fprintf(out, "%g\t%g\n", *data1++, *data2++);
}

void writepolar(FILE* out, const float* data, size_t len, double inc)
{  const float* data2 = data + len;
   // 1st line
   fprintf(out, "0\t%g\t0\n", *data++);
   len = 1;
   while (data < --data2)
      fprintf(out, "%g\t%g\t%g\n", len++*inc, *data++, *data2);
}

void writecomplex(FILE* out, const Complex* data, size_t len)
{  while (len--)
   {  fprintf(out, "%14g\t%14g\t%14g\t%14g\n", data->real(), data->imag(), abs(*data), arg(*data)*M_180_PI);
      ++data;
}  }

void write4complex(FILE* out, const Complex (* data)[4], size_t len)
{  while (len--)
   {  //Complex det = (*data)[0] * (*data)[3] - (*data)[1] * (*data)[2];
      fprintf(out, "%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\t%14g\n"
       , (*data)[0].real(), (*data)[0].imag(), (*data)[1].real(), (*data)[1].imag(), (*data)[2].real(), (*data)[2].imag(), (*data)[3].real(), (*data)[3].imag()
       , abs((*data)[0]), arg((*data)[0])*M_180_PI, abs((*data)[1]), arg((*data)[1])*M_180_PI, abs((*data)[2]), arg((*data)[2])*M_180_PI, abs((*data)[3]), arg((*data)[3])*M_180_PI);
       //, det.real(), det.imag());
      ++data;
}  }

void readcomplex(FILE*in, Complex * data, size_t len)
{  while (len--)
   {  double a,b;
      if (fscanf(in, "%lg%lg%*[^\n]", &a, &b) != 2)
         die("Failed to read complex data (%i).", errno);
      //(stderr, "%g\t%g\n", a,b);
      *data++ = Complex(a,b);
}  }

void read4complex(FILE*in, Complex (* data)[4], size_t len)
{  while (len--)
   {  double a,b,c,d,e,f,g,h;
      if (fscanf(in, "%lg%lg%lg%lg%lg%lg%lg%lg%*[^\n]", &a, &b, &c, &d, &e, &f, &g, &h) != 8)
         die("Failed to read complex data (%i).", errno);
      //(stderr, "%g\t%g\n", a,b);
      (*data)[0] = Complex(a,b);
      (*data)[1] = Complex(c,d);
      (*data)[2] = Complex(e,f);
      (*data)[3] = Complex(g,h);
      ++data;
}  }

/*static char* abbrev(const char* s, const char* token)
{  size_t l = strlen(s);
   return strnicmp(s, token, l) == 0 ? (char*)token + l : NULL;
}*/

static _cdecl void readint(const char* s, int* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%i%ln", r, &l) != 1 || l != strlen(s))
      die("Integer value expected");
}
static _cdecl void readintdef(const char* s, int* r, int d)
{  if (*s == 0)
      *r = d;
    else
   {  size_t l = (size_t)-1;
      if (sscanf(s, "%i%ln", r, &l) != 1 || l != strlen(s))
         die("Integer value expected");
}  }
static _cdecl void readfloat(const char* s, float* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%f%ln", r, &l) != 1 || l != strlen(s))
      die("Floating point value expected");
}
static _cdecl void readdouble(const char* s, double* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%lf%ln", r, &l) != 1 || l != strlen(s))
      die("Floating point value expected");
}

static _cdecl void readstring(const char* s, const char** cpp)
{  *cpp = s;
}

static _cdecl void setflag(const char* s, bool* r)
{  if (*s)
      die("Option does not have parameters");
   *r = true;
}

static _cdecl void setint(const char* s, int* r, int v)
{  if (*s)
      die("Option does not have parameters");
   *r = v;
}

static _cdecl void setbit(const char* s, int* r, int v)
{  if (*s)
      die("Option does not have parameters");
   *r |= v;
}

static _cdecl void readN(const char* s, int* r)
{  bool ex;
   if (ex = *s == '^')
      ++s;
   readint(s, r);
   if (ex)
      N = 1 << N;
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
   if (harmonic > 1)
   {  // compact to shorter FFT length Nh, Nh my be odd
      float* dp = outbuffer1; // real part of outbuffer1
      const float* sp = dp;
      const float* ep = sp + N/2;
      while ((sp += harmonic) < ep)
         *++dp = *sp;
      ep += N/2; // imaginary part of outbuffer1
      sp = ep - Nh/2 * harmonic;
      while (sp < ep)
      {  *++dp = *sp;
         sp += harmonic;
      }
      dp = outbuffer2; // real part of outbuffer2
      sp = dp;
      ep = sp + N/2;
      while ((sp += harmonic) < ep)
         *++dp = *sp;
      ep += N/2; // imaginary part of outbuffer2
      sp = ep - Nh/2 * harmonic;
      while (sp < ep)
      {  *++dp = *sp;
         sp += harmonic;
      }
   }
   // append constant zero to FFT result to simplify analysis
   outbuffer1[Nh] = 0;
   outbuffer2[Nh] = 0;
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

class FFTbin
{public:
   enum    StoreRet
   {  BelowMin,  // frequency less than fmin
      AboveMax,  // frequency above fmax
      Ready,     // calculated values available
      Aggregated // bin used for aggregation only
   };
 private:
   const double finc;

   int     binc;
   double  lphi; // last phase (for numerical derivative)

   double  af;
   Complex aU;
   Complex aI;
   Complex aZ;
   double  aD;
   double  cW;

 public:
   FFTbin(double finc) : finc(finc), binc(0), lphi(0) {}
   StoreRet StoreBin(int bin, Complex U, Complex I);

   double  f() const { return af; } // frequency
   Complex U() const { return aU; } // voltage, nominator or wanted signal
   Complex I() const { return aI; } // cuurrent, denominator or reference signal
   Complex Z() const { return aZ; } // impedance, quotient or relative signal
   double  D() const { return aD; } // group delay
   double  W() const { return cW; } // weight
};

FFTbin::StoreRet FFTbin::StoreBin(int bin, Complex U, Complex I)
{  // frequency
   double f = bin * finc;
   if (f < fmin)
      return BelowMin;
   if (f > fmax)
      return AboveMax;
   // calibration
   docal(bin, f, U, I);
   // phase correction
   I *= exp(Complex(0, linphase * f));
   // calc Y
   Complex Z(U/I);
   // group delay
   double D;
   {  double Zphi = arg(Z);
      D = (Zphi - lphi) / M_2PI;
      D -= floor(D+.5); // minimum phase
      //fprintf(stderr, "D: %f, %f\n", f, D);
      D /= finc;
      lphi = Zphi;
   }
   if (nophase)
      Z = abs(Z);
   // binsize
   if (fbinsc)
      binsz = (int)(bin * fbinsc + 1);
   if (binc == 0)
   {  // init
      af = f;
      aU = U;
      aI = I;
      aZ = Z;
      aD = D;
   } else
   {  af += f;
      aU += U;
      aI += I;
      aZ += Z;
      aD += D;
   }
   if (++binc != binsz)
      return Aggregated;
   af /= binc;
   aU /= binc;
   aI /= binc;
   aZ /= binc;
   aD /= binc;
   binc = 0;
   // weight
   cW = (*weightfn)(abs(aU), abs(aI), af);
   return Ready;
}

static void PrintFFTbin(FILE* dst, const FFTbin& calc)
{  fprintf(dst, "%12g %12g %12g %12g %12g " "%12g %12g %12g %12g %12g %12g\n",
   // f         |Hl|           phil                    |Hr|           phir
      calc.f(), abs(calc.U()), arg(calc.U())*M_180_PI, abs(calc.I()), arg(calc.I())*M_180_PI,
   // |Hl|/|Hr|      phil-phir               re               im
      abs(calc.Z()), arg(calc.Z())*M_180_PI, calc.Z().real(), calc.Z().imag(),
   // weight    delay
      calc.W(), calc.D());
}


static const struct ArgMap
{  char arg[8];
   void (_cdecl *func)(const char* rem, void* param, int iparam);
   void* param;
   int iparam;
} argmap[] = // must be sorted
{  {"ainc", (void(_cdecl*)(const char*,void*,int))&setflag, &incremental, true}
 , {"al",   (void(_cdecl*)(const char*,void*,int))&readint, &addloop}
 , {"bin",  (void(_cdecl*)(const char*,void*,int))&readint, &binsz}
 , {"ca" ,  (void(_cdecl*)(const char*,void*,int))&readintdef, &addch, 2}
 , {"df" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &datafile}
 , {"exec", (void(_cdecl*)(const char*,void*,int))&readstring, &execcmd}
 , {"famax",(void(_cdecl*)(const char*,void*,int))&readdouble, &famax}
 , {"famin",(void(_cdecl*)(const char*,void*,int))&readdouble, &famin}
 , {"fbin", (void(_cdecl*)(const char*,void*,int))&readdouble, &fbinsc}
 , {"fmax", (void(_cdecl*)(const char*,void*,int))&readdouble, &fmax}
 , {"fmin", (void(_cdecl*)(const char*,void*,int))&readdouble, &fmin}
 , {"fq" ,  (void(_cdecl*)(const char*,void*,int))&readdouble, &freq}
 , {"g2f" , (void(_cdecl*)(const char*,void*,int))&readstring, &gaindifffile}
 , {"gd" ,  (void(_cdecl*)(const char*,void*,int))&setint, &gainmode, 3}
 , {"gf" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &gainfile}
 , {"gg" ,  (void(_cdecl*)(const char*,void*,int))&setint, &gainmode, 2}
 , {"gr" ,  (void(_cdecl*)(const char*,void*,int))&setint, &gainmode, 1}
 , {"h/f",  (void(_cdecl*)(const char*,void*,int))&setint, &weightfn, (int)get1_fweight}
 , {"har",  (void(_cdecl*)(const char*,void*,int))&readint, &harmonic}
 , {"hd" ,  (void(_cdecl*)(const char*,void*,int))&setint, &weightfn, (int)getweightD}
 , {"he" ,  (void(_cdecl*)(const char*,void*,int))&setint, &weightfn, (int)getconstweight}
 , {"in" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &infile}
 , {"ln" ,  (void(_cdecl*)(const char*,void*,int))&readint, &loops, 1}
 , {"loop", (void(_cdecl*)(const char*,void*,int))&setint, &loops, INT_MAX}
 , {"lp" ,  (void(_cdecl*)(const char*,void*,int))&readint, &lpause}
 , {"lvl",  (void(_cdecl*)(const char*,void*,int))&readdouble, &noiselvl}
 , {"mfft", (void(_cdecl*)(const char*,void*,int))&setbit, &method, 1}
 , {"mpca", (void(_cdecl*)(const char*,void*,int))&setbit, &method, 2}
 , {"mxy",  (void(_cdecl*)(const char*,void*,int))&setbit, &method, 4}
 , {"n"  ,  (void(_cdecl*)(const char*,void*,int))&readN, &N}
 , {"pdc",  (void(_cdecl*)(const char*,void*,int))&readintdef, &purgech, 1}
 , {"phl",  (void(_cdecl*)(const char*,void*,int))&readdouble, &linphase}
 , {"phn",  (void(_cdecl*)(const char*,void*,int))&setint, &nophase, 1}
 , {"plot", (void(_cdecl*)(const char*,void*,int))&readstring, &plotcmd}
 , {"psa",  (void(_cdecl*)(const char*,void*,int))&readintdef, &discsamp, 1}
 , {"pte",  (void(_cdecl*)(const char*,void*,int))&readintdef, &disctrail, 1}
 , {"rf" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &rawfile}
 , {"rref", (void(_cdecl*)(const char*,void*,int))&readdouble, &rref}
 , {"scm",  (void(_cdecl*)(const char*,void*,int))&readint, &scalemode}
 , {"wd" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writedata, true}
 , {"wf" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &windowfile}
 , {"win",  (void(_cdecl*)(const char*,void*,int))&readintdef, &winfn, 2}
 , {"wr" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writeraw, true}
 , {"ww" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writewindow, true}
 , {"z2f" , (void(_cdecl*)(const char*,void*,int))&readstring, &zerodifffile}
 , {"zd" ,  (void(_cdecl*)(const char*,void*,int))&setint, &zeromode, 3}
 , {"zf" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &zerofile}
 , {"zg" ,  (void(_cdecl*)(const char*,void*,int))&setint, &zeromode, 2}
 , {"zn" ,  (void(_cdecl*)(const char*,void*,int))&setint, &normalize, 1}
 , {"zr" ,  (void(_cdecl*)(const char*,void*,int))&setint, &zeromode, 1}
};

static int searcharg(const char* arg, const char* elem)
{  return strnicmp(arg, elem, strlen(elem));
}

static void parsearg(const char* arg)
{  ArgMap* ap = (ArgMap*)bsearch(arg, argmap, sizeof argmap / sizeof *argmap, sizeof *argmap, (int (*)(const void*, const void*))&searcharg);
   if (ap == NULL)
      die("illegal option %s", arg);
   (*ap->func)(arg + strlen(ap->arg), ap->param, ap->iparam);
}

int main(int argc, char* argv[])
{  // parse cmdl
   while (--argc)
      parsearg(*++argv);

   if (N > N_MAX)
      die("FFT Length too large.");
   if (N * addch > N_MAX * CA_MAX)
      die("Input data length too large.");

   // prepare some global vars
   noiselvl_ = 1/noiselvl;
   freq /= addch;
   Nh = N / harmonic;
   linphase *= M_2PI;
   // allocate buffers
   if (method & 4)
   {  // reserver space for integrals and differentials too
      inbuffer1 = new float[max(N,3*Nh)];
      inbuffer2 = new float[max(N,3*Nh)];
      outbuffer1 = new float[max(N+1,3*Nh+1)];
      outbuffer2 = new float[max(N+1,3*Nh+1)];
   } else
   {  inbuffer1 = new float[N];
      inbuffer2 = new float[N];
      outbuffer1 = new float[N+1];
      outbuffer2 = new float[N+1];
   }
   wsums = new double[N_MAX/2+1];
   gain = new Complex[N_MAX/2+1];
   gainD = new Complex[N_MAX/2+1];
   zero = new Complex[N_MAX/2+1][4];
   zeroD = new Complex[N_MAX/2+1][4];
   // create plan
   //float in[N], tout[N], power_spectrum[N/2+1];
   P = rfftw_create_plan(N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
   PI = rfftw_create_plan(Nh, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);

   createwindow(window, winfn, N);
   // write window data
   FILE* tout;
   if (writewindow)
   {  tout = fopen(windowfile, "wt");
      if (tout == NULL)
         die("Failed to open %s.", windowfile);
      write1ch(tout, window, N);
      fclose(tout);
   }

   FILE* in;
   if (infile == NULL)
   {  // streaming
      in = stdin;
      _fsetmode(in, "b");
   } else
   {  in = fopen(infile, "rb");
      if (in == NULL)
         die("Failed to open input file.");
   }

   // discard first samples
   fread(inbuffertmp, sizeof *inbuffertmp, discsamp, in);

   memset(inbuffer1, 0, sizeof inbuffer1); // init woth 0 because of incremental mode
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
         die("Failed to open %s.", gainfile);
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
         die("Failed to open %s.", zerofile);
      read4complex(fz, zero, N/2);
      fclose(fz);
   }

   // operation loop
   int loop = loops;
   int addloops = 0;
   do
   {  if (fread(inbuffertmp, 2*sizeof *inbuffertmp * addch, N, in) != N)
         die("failed to read data");

      // reset min/max
      init();

      // write raw data
      if (writeraw)
      {  tout = fopen(rawfile, "wt");
         write2ch(tout, inbuffertmp, N * addch);
         fclose(tout);
      }

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
         die("Invalid scalmode");
      }
      // write raw status
      fprintf(stderr, "\nmin:\t%i\t%i\nmax:\t%i\t%i\n", minmax[0], minmax[2], minmax[1], minmax[3]);

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
         die("Invalid combination of measurement modes (%i), e.g. FFT and XY.", method);

       case 2: // PCA analysis
         {  PCA<5> pca;
            double data[6];
            float* U = inbuffer1 +1;
            float* I = inbuffer2 +1;
            data[2] = 1;   // offset
            data[3] = 0;   // integral
            data[4] = 0;   // linear
            data[5] = 0;   // differential
            int i = N -3;
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

            vectorscale(outbuffer1, sqrt(1./N)/addch, Nh);
            vectorscale(outbuffer2, sqrt(1./N)/addch, Nh);

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
            if (writedata)
            {  tout = fopen(datafile, "wt");
               if (tout == NULL)
                  die("Failed to create %s.", datafile);
            }
            // f l r l.arg l.ph r.arg r.ph
            const float* a1 = outbuffer1;
            const float* a2 = outbuffer2;
            const float* b1 = a1 + Nh;
            const float* b2 = a2 + Nh;

            FFTbin calc(freq/Nh);

            // 1st line
            for (size_t len = 0; a1 < b1; ++len, ++a1, ++a2, --b1, --b2)
            {  // do calculations and aggregations
               switch (calc.StoreBin(len, Complex(*a1, *b1), Complex(*a2, *b2)))
               {case FFTbin::BelowMin:
                case FFTbin::Aggregated:
                  continue;
                case FFTbin::AboveMax:
                  goto end1;
                default:; //case FFTbin::Ready:
               }
               // write
               if (writedata)
                  PrintFFTbin(tout, calc);

               if (calc.f() < famin && calc.f() >= famax)
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
          end1:
            if (writedata)
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

         }
         break;

       case 3: // FFT and then PCA
         {  // FFT
            dofft();

            vectorscale(outbuffer1, sqrt(1./N)/addch, Nh);
            vectorscale(outbuffer2, sqrt(1./N)/addch, Nh);

            // write data
            if (writedata)
            {  tout = fopen(datafile, "wt");
               if (tout == NULL)
                  die("Failed to create %s.", datafile);
            }
            // f l r l.arg l.ph r.arg r.ph
            const float* a1 = outbuffer1;
            const float* a2 = outbuffer2;
            const float* b1 = a1 + Nh;
            const float* b2 = a2 + Nh;

            PCA<2> pcaRe;
            PCA<3> pcaIm;
            double PCAdataRe[2];
            double PCAdataIm[3];
            // some values are const
            PCAdataRe[1] = 1;
            //PCAdataIm[1] = 1;

            FFTbin calc(freq/Nh);

            // 1st line
            for (size_t len = 0; a1 < b1; ++len, ++a1, ++a2, --b1, --b2)
            {  // do calculations and aggregations
               switch (calc.StoreBin(len, Complex(*a1, *b1), Complex(*a2, *b2)))
               {case FFTbin::BelowMin:
                case FFTbin::Aggregated:
                  continue;
                case FFTbin::AboveMax:
                  goto end3;
                default:; //case FFTbin::Ready:
               }
               // write
               if (writedata)
                  PrintFFTbin(tout, calc);

               if (calc.f() < famin && calc.f() >= famax)
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
          end3:
            if (writedata)
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
         }
         break;

       case 4: // XY mode
         {  // we need to do an FFT here, at least for the calibration
            dofft();
            // normalize (including retransformation)
            vectorscale(outbuffer1, sqrt(1./N/Nh)/addch, Nh);
            vectorscale(outbuffer2, sqrt(1./N/Nh)/addch, Nh);

            const double inc = freq/Nh;
            // U(f)
            float* a1 = outbuffer1;
            float* a2 = outbuffer2;
            float* b1 = a1 + Nh;
            float* b2 = a2 + Nh;
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
               a1[Nh] = C.real();
               b1[Nh] = C.imag();
               C = I/di;
               a2[Nh] = C.real();
               b2[Nh] = C.imag();
               // store differentials
               C = U*di;
               a1[2*Nh] = C.real();
               b1[2*Nh] = C.imag();
               C = I*di;
               a2[2*Nh] = C.real();
               b2[2*Nh] = C.imag();
            }
            // clear DC and nyquist frequencies of integral and diferential, since they do no allow 90¯ phase shift.
            outbuffer1[Nh] = outbuffer2[Nh] = 0;
            outbuffer1[2*Nh] = outbuffer2[2*Nh] = 0;
            outbuffer1[Nh+Nh/2] = outbuffer2[Nh+Nh/2] = 0;
            outbuffer1[2*Nh+Nh/2] = outbuffer2[2*Nh+Nh/2] = 0;

            // now do the inverse transform to get the corrected data back.
            rfftw_one(PI, outbuffer1, inbuffer1);
            rfftw_one(PI, outbuffer1 + Nh, inbuffer1 + Nh);
            rfftw_one(PI, outbuffer1 + 2*Nh, inbuffer1 + 2*Nh);
            rfftw_one(PI, outbuffer2, inbuffer2);
            rfftw_one(PI, outbuffer2 + Nh, inbuffer2 + Nh);
            rfftw_one(PI, outbuffer2 + 2*Nh, inbuffer2 + 2*Nh);

            // write data
            if (writedata)
            {  tout = fopen(datafile, "wt");
               if (tout == NULL)
                  die("Failed to create %s.", datafile);

               const float* Up = inbuffer1;
               const float* Ip = inbuffer2;
               for (int len = 0; len < Nh; ++len, ++Up, ++Ip)
                  fprintf(tout, "%12g %12g %12g %12g %12g %12g %12g\n",
                  // t         U    I    INT U  INT I  D U      D I
                     len/freq*harmonic, *Up, *Ip, Up[Nh], Ip[Nh], Up[2*Nh], Ip[2*Nh]);
               fclose(tout);
            }
         }
      }

      /*if (execcmd)
      {  system(execcmd);
      }*/
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
      double* wp = wsums + Nh/2+1;
      for (Complex* cp = gainD + Nh/2+1; --cp >= gainD; )
         *cp /= *--wp;  // scale 2 average
      FILE* fz;
      fz = fopen(gainmode == 3 ? gaindifffile : gainfile, "wb");
      if (fz == NULL)
         die("Failed to open %s.", gainmode == 3 ? gaindifffile : gainfile);
      writecomplex(fz, gainD, Nh/2+1);
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
            die("failed to read data");
         --loop;
      }
      puts("\nNow at part 2.");
      goto restart_zero;
    case 4: // part 2
    case 5:
      for (Complex (* cp)[4] = zeroD + Nh/2+1; --cp >= zeroD; )
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
         die("Failed to open %s.", zeromode == 5 ? zerodifffile : zerofile);
      write4complex(fz, zeroD, Nh/2);
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

