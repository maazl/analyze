#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <sys/builtin.h>
#include <stdarg.h>

#include <rfftw.h>
#include <complex>
using namespace std;
typedef complex<double> Complex;

#define M_2PI (2.*M_PI)
#define M_3PI (3.*M_PI)
#ifndef M_PI_2
#define M_PI_2 (M_PI*.5)
#endif
#define M_PI_180 (M_PI/180.)
#define M_180_PI (180./M_PI)

#define N_MAX (65536*8)
#define CA_MAX 2

// avoid name clash with math.h
#define fmin fmin__
#define fmax fmax__

// data buffers
static short inbuffertmp[2*CA_MAX*N_MAX];
static float inbuffer1[N_MAX];
static float inbuffer2[N_MAX];
static float outbuffer1[N_MAX+1];
static float outbuffer2[N_MAX+1];
static double sumbuffer1[N_MAX+1];
static double sumbuffer2[N_MAX+1];
static float window[N_MAX];
static Complex reference[N_MAX+2];


void die(const char* msg, ...)
{  va_list va;
   va_start(va, msg);
   vprintf(msg, va);
   va_end(va);
   putchar('\n');
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

static const double minval = 1E-20;

static void vectorscale(float* data, double factor, size_t len)
{  while (len--)
      *data++ *= factor;
}

static void vectorscalecopy(float* dst, const double* src, double factor, size_t len)
{  while (len--)
      *dst++ = *src++ * factor;
}

static void vectoradd(float* dst, float* src, size_t len)
{  while (len--)
      *dst++ += *src++;
}

static void vectoradd(double* dst, float* src, size_t len)
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
static bool       swapbytes   = false;
static int        winfn       = 0;
static double     freq        = 44100;
static double     fmin        = -1;
static double     fmax        = 1E99;
static double     famin       = 1;
static double     famax       = 1E99;
static bool       writeraw    = false;
static bool       writedata   = false;
static bool       writewindow = false;
static int        purgech     = 0;
static int        discsamp    = 0;
static int        addch       = 1;
static int        loops       = 1;
static int        binsz       = 1;
static double     fbinsc      = 0;
static bool       nophaseA    = false;
static bool       nophaseB    = false;
static int        refmode     = 0;
static const char* infile     = NULL;
static const char* execcmd    = NULL;
static const char* plotcmd    = NULL;


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
      {  d1 += storeminmax(src[0], minmax);
         d2 += storeminmax(src[1], minmax+2);
         src += 2;
      } while (--i);
      *dst1++ = d1 * gainadj[0];
      *dst2++ = d2 * gainadj[1];
      --len;
   }
}

static void shortx2float2(float* dst1, float* dst2, const short* src, size_t len)
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

static void shortx2float2window(float* dst1, float* dst2, const short* src, const float* win, size_t len)
{  while (len)
   {  double d1 = 0;
      double d2 = 0;
      int i = addch;
      do
      {  d1 += storeminmax(src[0], minmax);
         d2 += storeminmax(src[1], minmax+2);
         src += 2;
      } while (--i);
      *dst1++ = d1 * *win * gainadj[0];
      *dst2++ = d2 * *win++ * gainadj[1];
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

static void polar2complex(float* data, size_t len)
{  float *data2 = data + len;
   while (++data < --data2)
   {  register double r = *data;
      *data = r * cos(*data2);
      *data2 = r * sin(*data2);
   }
}

static void init()
{  minmax[0] = INT_MAX;
   minmax[1] = INT_MIN;
   minmax[2] = INT_MAX;
   minmax[3] = INT_MIN;
}

static void outdata(const float* left, const float* right)
{
   const double inc = freq/N;
   // f l r l.arg l.ph r.argr.ph
   const float* a1 = left;
   const float* a2 = right;
   const float* b1 = a1 + N;
   const float* b2 = a2 + N;
   const Complex* rp = reference;

   double mf;
   Complex m1, m2;

   FILE* tout = fopen("data.dat", "wt");
   if (tout == NULL)
      die("Failed to create data.dat.");

   // 1st line
   int binc = 0;
   for (size_t len = 0; a1 < b1; ++len, ++a1, ++a2, --b1, --b2, rp += 2)
   {  // prepare
      double f = len*inc;
      if (f <= fmin)
         continue;
      if (f > fmax)
         break;
      Complex c1 = Complex(*a1, *b1);
      Complex c2 = Complex(*a2, *b2);
      // reference
      if (refmode == 1)
      {  c1 /= rp[0];
         c2 /= rp[1];
      }
      // binsize
      if (fbinsc)
         binsz = (int)(f * fbinsc + 1);
      if (nophaseB)
      {  if (binc == 0)
         {  // init
            mf = f;
            m1 = Complex(abs(c1), arg(c1)); // pfusch
            m2 = Complex(abs(c2), arg(c2));
         } else
         {  mf += f;
            m1 += Complex(abs(c1), arg(c1));
            m2 += Complex(abs(c2), arg(c2));
         }
         if (++binc != binsz)
            continue;
         mf /= binc;
         m1 /= binc;
         m2 /= binc;
         m1 = Complex(m1.real() * cos(m1.imag()), m1.real() * sin(m1.imag()));
         m2 = Complex(m2.real() * cos(m2.imag()), m2.real() * sin(m2.imag()));
         binc = 0;
      } else
      {  if (binc == 0)
         {  // init
            mf = f;
            m1 = c1;
            m2 = c2;
         } else
         {  mf += f;
            m1 += c1;
            m2 += c2;
         }
         if (++binc != binsz)
            continue;
         mf /= binc;
         m1 /= binc;
         m2 /= binc;
         binc = 0;
      }
      // prepare
      Complex Z = m2/m1;
      // write
      fprintf(tout, "%12g%12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
      // f   |l|      phi l             |Hr|     phi r
         mf, abs(m1), arg(m1)*M_180_PI, abs(m2), arg(m2)*M_180_PI,
      // |r/l|   pli r/l          la         lb         ra         rb         Za        Zb
         abs(Z), arg(Z)*M_180_PI, m1.real(), m1.imag(), m2.real(), m2.imag(), Z.real(), Z.imag());

   }
   fclose(tout);
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

static void readcomplex(FILE*in, Complex* data, size_t len)
{  while (len--)
   {  double a,b;
      if (fscanf(in, "%lg%lg%*[^\n]", &a, &b) != 2)
         die("Failed to read complex data (%i).", errno);
      //(stderr, "%g\t%g\n", a,b);
      *data++ = Complex(a,b);
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

static const struct ArgMap
{  char arg[8];
   void (_cdecl *func)(const char* rem, void* param, int iparam);
   void* param;
   int iparam;
} argmap[] = // must be sorted
{  {"ap" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &nophaseA, true}
 , {"bin",  (void(_cdecl*)(const char*,void*,int))&readint, &binsz, 0}
 , {"bp" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &nophaseB, true}
 , {"ca" ,  (void(_cdecl*)(const char*,void*,int))&readintdef, &addch, 2}
 , {"exec", (void(_cdecl*)(const char*,void*,int))&readstring, &execcmd, 0}
 , {"famax",(void(_cdecl*)(const char*,void*,int))&readdouble, &famax, 0}
 , {"famin",(void(_cdecl*)(const char*,void*,int))&readdouble, &famin, 0}
 , {"fbin", (void(_cdecl*)(const char*,void*,int))&readdouble, &fbinsc, 0}
 , {"fmax", (void(_cdecl*)(const char*,void*,int))&readdouble, &fmax, 0}
 , {"fmin", (void(_cdecl*)(const char*,void*,int))&readdouble, &fmin, 0}
 , {"fq" ,  (void(_cdecl*)(const char*,void*,int))&readdouble, &freq, 0}
 , {"in" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &infile, 0}
 , {"ln" ,  (void(_cdecl*)(const char*,void*,int))&readint, &loops, 1}
 , {"loop", (void(_cdecl*)(const char*,void*,int))&setint, &loops, INT_MAX}
 , {"n"  ,  (void(_cdecl*)(const char*,void*,int))&readN, &N, 0}
 , {"pdc",  (void(_cdecl*)(const char*,void*,int))&readintdef, &purgech, 1}
 , {"plot", (void(_cdecl*)(const char*,void*,int))&readstring, &plotcmd, 0}
 , {"psa",  (void(_cdecl*)(const char*,void*,int))&readintdef, &discsamp, 1}
 , {"rg" ,  (void(_cdecl*)(const char*,void*,int))&setint, &refmode, 2}
 , {"rr" ,  (void(_cdecl*)(const char*,void*,int))&setint, &refmode, 1}
 , {"wd" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writedata, true}
 , {"win",  (void(_cdecl*)(const char*,void*,int))&readintdef, &winfn, 2}
 , {"wr" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writeraw, true}
 , {"ww" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writewindow, true}
 , {"xb" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &swapbytes, true}
};

static int searcharg(const char* arg, const char* elem)
{  return strnicmp(arg, elem, strlen(elem));
}

static void parsearg(const char* arg)
{  ArgMap* ap = (ArgMap*)bsearch(arg, argmap, sizeof argmap / sizeof *argmap, sizeof *argmap, (int (*)(const void*, const void*))&searcharg);
   if (ap == NULL)
      die("illegal option");
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

   freq /= addch;
   // create plan
   //float in[N], tout[N], power_spectrum[N/2+1];
   rfftw_plan p = rfftw_create_plan(N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
   // append constant zero to FFT result to simplify analysis
   outbuffer1[N] = 0;
   outbuffer2[N] = 0;

   createwindow(window, winfn, N);
   // write window data
   FILE* tout;
   if (writewindow)
   {  tout = fopen("window.dat", "wt");
      if (tout == NULL)
         die("Failed to open window.dat");
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

   memset(sumbuffer1, 0, sizeof sumbuffer1);
   memset(sumbuffer2, 0, sizeof sumbuffer2);

   // prepare refmode
   switch (refmode)
   {case 1: // read
      FILE* fz;
      fz = fopen("ref.dat", "r");
      if (fz == NULL)
         die("Failed to open ref.dat.");
      readcomplex(fz, reference, N+2);
      fclose(fz);
      break;
    case 2: // generate
      memset(reference, 0, sizeof reference);
   }

   // operation loop
   int loop = loops;
   do
   {  if (fread(inbuffertmp, 2*sizeof *inbuffertmp * addch, N, in) != N)
         die("failed to read data");

      init();

      // write raw data
      if (writeraw)
      {  tout = fopen("raw.dat", "wt");
         write2ch(tout, inbuffertmp, N * addch);
         fclose(tout);
      }

      if (winfn)
      {  if (swapbytes)
            shortx2float2window(inbuffer1, inbuffer2, inbuffertmp, window, N);
          else
            short2float2window(inbuffer1, inbuffer2, inbuffertmp, window, N);
      } else
      {  if (swapbytes)
            shortx2float2(inbuffer1, inbuffer2, inbuffertmp, N);
          else
            short2float2(inbuffer1, inbuffer2, inbuffertmp, N);
      }

      // write raw status
      fprintf(stderr, "\nmin:\t%i\t%i\nmax:\t%i\t%i\n", minmax[0], minmax[2], minmax[1], minmax[3]);

      // write raw data
      /*if (writeraw)
      {  tout = fopen("raw.dat", "wt");
         write2ch(tout, inbuffer1, inbuffer2, N);
         fclose(tout);
      }*/

      // FFT
      rfftw_one(p, inbuffer1, outbuffer1);
      rfftw_one(p, inbuffer2, outbuffer2);

      // purge DC (senseless without DC-coupling)
      if (purgech)
      {  static const double minscale = 1E-15;
         outbuffer1[0] *= minscale;
         outbuffer2[0] *= minscale;
         for (int i = purgech; --i;)
         {  outbuffer1[i] *= minscale;
            outbuffer2[i] *= minscale;
            outbuffer1[N-i] *= minscale;
            outbuffer2[N-i] *= minscale;
      }  }
      vectorscale(outbuffer1, sqrt(1./N)/addch, N);
      vectorscale(outbuffer2, sqrt(1./N)/addch, N);

      if (writedata)
         outdata(outbuffer1, outbuffer2);

      if (plotcmd)
      {  // for gnuplot!
         puts(plotcmd);
         fflush(stdout);
      }

      // calc sum
      if (nophaseA)
      {  complex2polar(outbuffer1, N);
         complex2polar(outbuffer2, N);
      }
      vectoradd(sumbuffer1, outbuffer1, N);
      vectoradd(sumbuffer2, outbuffer2, N);

   } while (--loop);
   // close stdin to signal playrec
   fclose(in);

   vectorscalecopy(outbuffer1, sumbuffer1, 1./loops, N);
   vectorscalecopy(outbuffer2, sumbuffer2, 1./loops, N);
   if (nophaseA)
   {  polar2complex(outbuffer1, N);
      polar2complex(outbuffer2, N);
   }

   outdata(outbuffer1, outbuffer2);

   if (refmode == 2) // write
   {  Complex* rp = reference;
      const float* a1 = outbuffer1;
      const float* a2 = outbuffer2;
      const float* b1 = a1 + N;
      const float* b2 = a2 + N;
      for (int i = N; --i;)
      {  *rp++ = Complex(*a1++, *b1--);
         *rp++ = Complex(*a2++, *b2--);
      }

      FILE* fz;
      fz = fopen("ref.dat", "wb");
      if (fz == NULL)
         die("Failed to open ref.dat.");
      writecomplex(fz, reference, N+2);
      fclose(fz);
   }

   return 0;
}
