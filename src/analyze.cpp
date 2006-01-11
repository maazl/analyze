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
static float      gainadj[2]  = {1,1};
static int        N           = 8192;
static double     noiselvl    = 1;
static int        winfn       = 0;
static double     freq        = 44100;
static double     fmin        = -1;
static double     fmax        = 1E99;
static double     famin       = 1;
static double     famax       = 1E99;
static bool       writeraw    = false;
static bool       writedata   = false;
static bool       writewindow = false;
static int        method      = 0;
static int        purgech     = 0;
static int        discsamp    = 0;
static bool       disctrail   = false;
static int        addch       = 1;
static double     rref        = 1;
static int        scalemode   = 1;
static int        loops       = 1;
static int        lpause      = 10;
static int        zeromode    = 0; // 0 = none, 1 = read, 2 = generate, 3 = generatedelta, 4 = generate part 2, 5 = generatedelta part 2
static int        gainmode    = 0; // 0 = none, 1 = read, 2 = generate, 3 = generatedelta
static bool       normalize   = false;
static int        binsz       = 1;
static double     fbinsc      = 0;
static int        harmonic    = 1;
static const char* datafile   = "data.dat";
static const char* zerofile   = "zero.dat";
static const char* zerodifffile = "zeroD.dat";
static const char* gainfile   = "gain.dat";
static const char* gaindifffile = "gainD.dat";
static const char* rawfile    = "raw.dat";
static const char* windowfile = "window.dat";
static const char* infile     = NULL;
static const char* execcmd    = NULL;
static const char* plotcmd    = NULL;
static double     (*weightfn)(double, double, double) = getweight;
// internal vars
static int        Nh;   // FFT length after harmonic extraction
static rfftw_plan P;
static rfftw_plan PI;


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

static void short2float2window(float* dst1, float* dst2, const short* src, const float* win, size_t len)
{  while (len)
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

static void short2floatDwindow(float* dst1, float* dst2, const short* src, const float* win, size_t len)
{  while (len)
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

/* Calibration */
static double docal(int len, double f, Complex& U, Complex& I)
{  /* introduce error */
   /*const Complex c11(1.1, .1);
   const Complex c12(.1, .01);
   const Complex c21(.3, .02);
   const Complex c22(.95, .2);
   const Complex c11(1, 0);
   const Complex c12(0, 0);
   const Complex c21(0, 0);
   const Complex c22(1, 0);
   {  Complex t = U;
      U = c11 * U + c12 * I;
      I = c21 * t + c22 * I;
   }*/
   double weight;
   switch (gainmode)
   {case 1: // read
      U /= gain[len];
      /*t = I * gainD[len];
      I -= U * gainD[len];
      U -= t;*/
      break;
    case 2: // write
    case 3:
      wsums[len] += weight = sqrt((*weightfn)(abs(U), abs(I), f));
      U /= gain[len];
      gainD[len] += I/U * weight;
   }
   if (normalize)
   {  Complex s = 1./(U + I);
      U *= s;
      I *= s;
   }
   if (zeromode & 1)
   {  Complex* cp = zero[len];
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
      zeroD[len][0] += U;
      zeroD[len][2] += I;
      break;
    case 4: // write part 2
    case 5:
      zeroD[len][1] += U;
      zeroD[len][3] += I;
   }
   return weight;
}


static const struct ArgMap
{  char arg[8];
   void (_cdecl *func)(const char* rem, void* param, int iparam);
   void* param;
   int iparam;
} argmap[] = // must be sorted
{  {"bin",  (void(_cdecl*)(const char*,void*,int))&readint, &binsz, 0}
 , {"ca" ,  (void(_cdecl*)(const char*,void*,int))&readintdef, &addch, 2}
 , {"df" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &datafile, 0}
 , {"exec", (void(_cdecl*)(const char*,void*,int))&readstring, &execcmd, 0}
 , {"famax",(void(_cdecl*)(const char*,void*,int))&readdouble, &famax, 0}
 , {"famin",(void(_cdecl*)(const char*,void*,int))&readdouble, &famin, 0}
 , {"fbin", (void(_cdecl*)(const char*,void*,int))&readdouble, &fbinsc, 0}
 , {"fmax", (void(_cdecl*)(const char*,void*,int))&readdouble, &fmax, 0}
 , {"fmin", (void(_cdecl*)(const char*,void*,int))&readdouble, &fmin, 0}
 , {"fq" ,  (void(_cdecl*)(const char*,void*,int))&readdouble, &freq, 0}
 , {"g2f" , (void(_cdecl*)(const char*,void*,int))&readstring, &gaindifffile, 0}
 , {"gd" ,  (void(_cdecl*)(const char*,void*,int))&setint, &gainmode, 3}
 , {"gf" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &gainfile, 0}
 , {"gg" ,  (void(_cdecl*)(const char*,void*,int))&setint, &gainmode, 2}
 , {"gr" ,  (void(_cdecl*)(const char*,void*,int))&setint, &gainmode, 1}
 , {"h/f",  (void(_cdecl*)(const char*,void*,int))&setint, &weightfn, (int)get1_fweight}
 , {"har",  (void(_cdecl*)(const char*,void*,int))&readint, &harmonic, 0}
 , {"hd" ,  (void(_cdecl*)(const char*,void*,int))&setint, &weightfn, (int)getweightD}
 , {"he" ,  (void(_cdecl*)(const char*,void*,int))&setint, &weightfn, (int)getconstweight}
 , {"in" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &infile, 0}
 , {"ln" ,  (void(_cdecl*)(const char*,void*,int))&readint, &loops, 1}
 , {"loop", (void(_cdecl*)(const char*,void*,int))&setint, &loops, INT_MAX}
 , {"lp" ,  (void(_cdecl*)(const char*,void*,int))&readint, &lpause, 0}
 , {"lvl",  (void(_cdecl*)(const char*,void*,int))&readdouble, &noiselvl, 0}
 , {"mfft", (void(_cdecl*)(const char*,void*,int))&setbit, &method, 1}
 , {"mpca", (void(_cdecl*)(const char*,void*,int))&setbit, &method, 2}
 , {"mxy",  (void(_cdecl*)(const char*,void*,int))&setbit, &method, 4}
 , {"n"  ,  (void(_cdecl*)(const char*,void*,int))&readN, &N, 0}
 , {"pdc",  (void(_cdecl*)(const char*,void*,int))&readintdef, &purgech, 1}
 , {"plot", (void(_cdecl*)(const char*,void*,int))&readstring, &plotcmd, 0}
 , {"psa",  (void(_cdecl*)(const char*,void*,int))&readintdef, &discsamp, 1}
 , {"pte",  (void(_cdecl*)(const char*,void*,int))&readintdef, &disctrail, 1}
 , {"rf" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &rawfile, 0}
 , {"rref", (void(_cdecl*)(const char*,void*,int))&readdouble, &rref, 0}
 , {"scm",  (void(_cdecl*)(const char*,void*,int))&readint, &scalemode, 0}
 , {"wd" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writedata, true}
 , {"wf" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &windowfile, 0}
 , {"win",  (void(_cdecl*)(const char*,void*,int))&readintdef, &winfn, 2}
 , {"wr" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writeraw, true}
 , {"ww" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writewindow, true}
 , {"z2f" , (void(_cdecl*)(const char*,void*,int))&readstring, &zerodifffile, 0}
 , {"zd" ,  (void(_cdecl*)(const char*,void*,int))&setint, &zeromode, 3}
 , {"zf" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &zerofile, 0}
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
   do
   {  if (fread(inbuffertmp, 2*sizeof *inbuffertmp * addch, N, in) != N)
         die("failed to read data");

      init();

      // write raw data
      if (writeraw)
      {  tout = fopen(rawfile, "wt");
         write2ch(tout, inbuffertmp, N * addch);
         fclose(tout);
      }

      //
      switch (scalemode)
      {case 1:
         if (winfn)
            short2float2window(inbuffer1, inbuffer2, inbuffertmp, window, N);
          else
            short2float2(inbuffer1, inbuffer2, inbuffertmp, N);
         break;
       case 2:
         if (winfn)
            short2floatDwindow(inbuffer1, inbuffer2, inbuffertmp, window, N);
          else
            short2floatD(inbuffer1, inbuffer2, inbuffertmp, N);
         break;
       case 3:
         if (winfn)
            short2float2window(inbuffer2, inbuffer1, inbuffertmp, window, N);
          else
            short2float2(inbuffer2, inbuffer1, inbuffertmp, N);
         break;
       default:
         die("Invalid scalmode");
      }
      // write raw status
      fprintf(stderr, "\nmin:\t%i\t%i\nmax:\t%i\t%i\n", minmax[0], minmax[2], minmax[1], minmax[3]);

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
               data[5] = I[-1] + I[0] - I[1] - I[2];
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

            const double inc = freq/Nh;
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

            double af;
            Complex aU;
            Complex aI;
            Complex aZ;

            // 1st line
            int binc = 0;
            for (size_t len = 0; a1 < b1; ++len, ++a1, ++a2, --b1, --b2)
            {  // calc
               double f = len*inc;
               if (f < fmin)
                  continue;
               if (f > fmax)
                  break;
               Complex U(*a1, *b1);
               Complex I(*a2, *b2);
               // calibration
               double weight = docal(len, f, U, I);
               // calc Y
               Complex Z(U/I);
               // binsize
               {  if (fbinsc) // dynamic binsize
                     binsz = (int)(f * fbinsc + 1);
                  if (binc == 0)
                  {  // init
                     af = f;
                     aU = U;
                     aI = I;
                     aZ = Z;
                  } else
                  {  af += f;
                     aU += U;
                     aI += I;
                     aZ += Z;
                  }
                  if (++binc != binsz)
                     continue;
                  af /= binc;
                  aU /= binc;
                  aI /= binc;
                  aZ /= binc;
                  binc = 0;
               }
               // prepare
               double absU = abs(aU);
               double absI = abs(aI);
               weight = (*weightfn)(absU, absI, af);
               // write
               if (writedata)
                  fprintf(tout, "%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
                  // f  |Hl|  phil             |Hr|  phir
                     af, absU, arg(aU)*M_180_PI, absI, arg(aI)*M_180_PI,
                     //f, *a1 , *b1            , *a2 , *b2            ,
                  // |Hl|/|Hr|  phil-phir        re        im
                     abs(aZ), arg(aZ)*M_180_PI, aZ.real(), aZ.imag(), weight);

               if (af < famin && af >= famax)
                  continue;
               // average
               ++nsum;
               wsum += weight;
               // resistivity
               Rsum += weight * aZ.real();
               R2sum += weight * sqr(aZ.real());
               // L & C
               L2sum += weight * sqr(af);
               // LCsum += weight; == wsum
               C2sum += weight / sqr(af);
               Lsum -= weight * af * aZ.imag();
               Csum -= weight / af * aZ.imag();
               d2sum = weight * sqr(aZ.imag());
            }
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

            const double inc = freq/Nh;
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

            double af;
            Complex aU;
            Complex aI;
            Complex aZ;

            // 1st line
            int binc = 0;
            for (size_t len = 0; a1 < b1; ++len, ++a1, ++a2, --b1, --b2)
            {  // calc
               double f = len*inc;
               if (f < fmin)
                  continue;
               if (f > fmax)
                  break;
               Complex U(*a1, *b1);
               Complex I(*a2, *b2);
               // calibration
               double weight = docal(len, f, U, I);
               // calc Y
               Complex Z(U/I);
               // binsize
               {  if (fbinsc)
                     binsz = (int)(f * fbinsc + 1);
                  if (binc == 0)
                  {  // init
                     af = f;
                     aU = U;
                     aI = I;
                     aZ = Z;
                  } else
                  {  af += f;
                     aU += U;
                     aI += I;
                     aZ += Z;
                  }
                  if (++binc != binsz)
                     continue;
                  af /= binc;
                  aU /= binc;
                  aI /= binc;
                  aZ /= binc;
                  binc = 0;
               }
               // prepare
               double absU = abs(aU);
               double absI = abs(aI);
               weight = (*weightfn)(absU, absI, af);
               // write
               if (writedata)
                  fprintf(tout, "%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
                  // f   |Hl|  phil              |Hr|  phir
                     af, absU, arg(aU)*M_180_PI, absI, arg(aI)*M_180_PI,
                  // |Hl|/|Hr| phil-phir        re         im
                     abs(aZ), arg(aZ)*M_180_PI, aZ.real(), aZ.imag(), weight);

               if (af < famin && af >= famax)
                  continue;
               // component analysis
               PCAdataRe[0] = aZ.real();
               //PCAdataRe[2] = 1/af;
               //PCAdataRe[3] = f;
               PCAdataIm[0] = aZ.imag(); // // fit imaginary part in conductivity
               PCAdataIm[1] = 1/af;
               PCAdataIm[2] = af;
               //PCAdataIm[3] = 1/af;
               //printf("Re: %12g %12g %12g %12g\n", PCAdataRe[0], PCAdataRe[1], PCAdataRe[2], weight);

               // add values
               pcaRe.Store(PCAdataRe, weight);
               pcaIm.Store(PCAdataIm, weight);

            }
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

