#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <386/builtin.h>
#include <stdarg.h>

#ifdef __OS2__
#define INCL_BASE
#include <os2.h>
#endif

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
#define M_sqrtPI (1.772453851)

#define N_MAX (65536*8)
#define CA_MAX 2
#define SYNC_FIR 64


int __gxx_personality_v0; // gcc @õ$%&!


// data buffers
static short refbuffer[N_MAX*2];
static short inprebuffer[CA_MAX*(N_MAX+SYNC_FIR)];
static short* const inbuffer = inprebuffer + CA_MAX*SYNC_FIR;
static float outprebuffer[N_MAX+2*SYNC_FIR];
static float* const outbuffer = outprebuffer + 2*SYNC_FIR;

static bool termrq = false;
static FILE* resh = NULL;

void die(const char* msg, ...)
{  termrq = true;
   va_list va;
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

static const double minval = 1E-20;


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
static int        N           = 8192;
static double     freq        = 48000;
static double     fmin        = 20;
static double     fmax        = 20000;
static double     fstep       = 1.05946309;
static int        discardsamp = 0;
static int        syncsamp    = 50000;
static int        overlap     = 1000;
static double     synclevel   = 10000;
static double     syncphase   = 1;
static bool       writeraw    = false;
static bool       writedata   = false;
static int        loops       = 1;
static int        mode        = 3; // 1 = generate, 2 = analyze 3 = both
static int        zeromode    = 0; // 0 = none, 1 = read, 2 = generate
static int        gainmode    = 0; // 0 = none, 1 = read, 2 = generate, 3 = generatedelta
static bool       verbose     = false;
static const char* infile     = NULL;
static const char* execcmd    = NULL;

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

static void short2float2(float* dst1, float* dst2, const short* src, size_t len)
{  while (len)
   {  *dst1++ = storeminmax(fromraw(src[0]), minmax);
      *dst2++ = storeminmax(fromraw(src[1]), minmax+2);
      --len;
      src += 2;
   }
}

static inline double abs(double d1, double d2)
{  return sqrt(sqr(d1) + sqr(d2));
}

void init()
{  minmax[0] = INT_MAX;
   minmax[1] = INT_MIN;
   minmax[2] = INT_MAX;
   minmax[3] = INT_MIN;
}

/*void write1ch(FILE* out, const float* data, size_t len)
{  while (len--)
      fprintf(out, "%g\n", *data++);
}*/

void write2ch(FILE* out, const short* data, size_t len)
{  while (len--)
   {  fprintf(out, "%i\t%i\n", fromraw(data[0]), fromraw(data[1]));
      data += 2;
   }
}

/*void write2ch(FILE* out, const float* data, size_t len)
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

void readcomplex(FILE*in, Complex* data, size_t len)
{  while (len--)
   {  double a,b;
      if (fscanf(in, "%lg%lg%*[^\n]", &a, &b) != 2)
         die("Failed to read complex data (%i).", errno);
      //(stderr, "%g\t%g\n", a,b);
      *data++ = Complex(a,b);
}  }*/

static void fwriteexact(const void* buffer, size_t size, size_t n, FILE* f)
{  while (n)
   {  size_t r = fwrite(buffer, size, n, f);
      if (r <= 0)
         die("Failed to write %lu blocks a %lu bytes", n, size);
      n -= r;
      (const char*&)buffer += r;
}  }

static void freadexact(void* buffer, size_t size, size_t n, FILE* f)
{  while (n)
   {  size_t r = fread(buffer, size, n, f);
      if (r <= 0)
         die("Failed to read %lu blocks a %lu bytes", n, size);
      //fwrite(buffer, size, r, stdout);
      n -= r;
      (char*&)buffer += r;
   }
}



// synchronize
static void syncinit()
{  memset(inprebuffer, 0, sizeof inprebuffer);
   memset(outprebuffer, 0, sizeof outprebuffer);
}

static bool issync(float* dp)
{
   static int synccount = 0;
   if (verbose)
   {  static int cnt = 0;
      fprintf(stderr, "# Res: %2x %- 7f %- 7f %i %i %i %i %- 7f %- 7f %- 7f\n",
       ++cnt % SYNC_FIR, dp[0], M_180_PI*dp[1],
       dp[0] >= synclevel, dp[-2*SYNC_FIR] >= synclevel, M_180_PI*abs(fmod(dp[1] - dp[1-2*SYNC_FIR] + M_2PI, M_2PI) - M_PI) <= syncphase, synccount,
       dp[-2*SYNC_FIR], M_180_PI*abs(fmod(dp[1] - dp[1-2*SYNC_FIR] + M_2PI, M_2PI) - M_PI), M_180_PI*abs(fmod(dp[1] - dp[1-2*SYNC_FIR] + M_2PI, M_2PI) - M_PI));
   }

   if ( dp[0] >= synclevel
       && dp[-2*SYNC_FIR] >= synclevel
       && M_180_PI*abs(fmod(dp[1] - dp[1-2*SYNC_FIR] + M_2PI, M_2PI) - M_PI) <= syncphase )
      return ++synccount == SYNC_FIR>>2;
   synccount = 0;
   return false;
}

static int synchronize(int len)
{  // FIR Filter
   static int rem = 0;
   const short* sp = inbuffer;
   float* dp = outbuffer;
   const short* se = sp + 2*len;
   sp += rem;
   while (sp < se)
   {  const short* sp2 = sp;
      int a = 0;
      int b = 0;
      for (int l = SYNC_FIR>>2; l; l--)
      {  a += fromraw(sp2[0]) - fromraw(sp2[-4]);
         b += fromraw(sp2[-2]) - fromraw(sp2[-6]);
         sp2 -= 8;
      }
      dp[0] = abs(a, b);
      dp[1] = atan2((float)b, a);
      //fprintf(stderr, "# Data: %4.4x %4.4x %4.4x %4.4x\t", fromraw(sp[0]), fromraw(sp[-2]), fromraw(sp[-4]), fromraw(sp[-6]));
      if (issync(dp))
         return (sp - inbuffer) >> 1;
      dp += 2;
      sp += 8;
   }
   rem = sp - se;
   // save history
   memcpy(inprebuffer, sp-2*SYNC_FIR, 2*SYNC_FIR * sizeof(short));
   memcpy(outprebuffer, dp-2*SYNC_FIR, 2*SYNC_FIR * sizeof(float));
   return -1;
}


// analysis
static double ana[10];

static void analyze(int fi)
{  init();
   memset(ana, 0, sizeof ana);
   const double fs = M_2PI*fi/N;
   const short* sp = inbuffer;
   for (int i = 0; i < N; ++i)
   {  double s = sin(fs*i);
      double c = cos(fs*i);
      // calculate sin/cos sums
      double v = storeminmax(fromraw(sp[0]), minmax);
      ana[0] += c * v;
      ana[1] += s * v;
      ana[6] += v;
      ana[8] += sqr(v);
      v = storeminmax(fromraw(sp[1]), minmax+2);
      ana[2] += c * v;
      ana[3] += s * v;
      ana[7] += v;
      ana[9] += sqr(v);
      sp += 2;
   }
   ana[0] /= N/M_SQRT2;
   ana[1] /= N/M_SQRT2;
   ana[2] /= N/M_SQRT2;
   ana[3] /= N/M_SQRT2;
   ana[6] /= N;
   ana[7] /= N;
   ana[8] /= N;
   ana[9] /= N;
   double s = sqr(ana[0]) + sqr(ana[1]);
   ana[4] = (ana[0]*ana[2] + ana[1]*ana[3]) / s;
   ana[5] = (ana[0]*ana[3] - ana[1]*ana[2]) / s;
   //printf("X: %f\t%f\t%f\t%f\n", ana[6], ana[7], sqrt(ana[8]), sqrt(ana[9]));
   ana[6] = sqrt(ana[8] - s - sqr(ana[6]));
   ana[7] = sqrt(ana[9] - sqr(ana[2]) - sqr(ana[3]) - sqr(ana[7]));
}

static void doanalysis(void*)
{  FILE* in;
   if (infile == NULL)
   {  // streaming
      in = stdin;
      _fsetmode(in, "b");
   } else
   {  in = fopen(infile, "rb");
      if (in == NULL)
         die("Failed to open input file.");
   }

   //_fsetmode(stdout, "b"); // debug

   // skip initial samples
   freadexact(inbuffer, 2*sizeof(short), discardsamp, in);

   // synchronize
   int i;
   int n = 0;
   for(;;)
   {  if (n >= syncsamp)
         die("Failed to syncronize.");
      freadexact(inbuffer, 2*sizeof(short), syncsamp/4, in);
      /*FILE* fs = fopen("sync.dat", "a");
      write2ch(fs, inbuffer, syncsamp/4);
      fclose(fs);*/
      i = synchronize(syncsamp>>2);
      if (i >= 0)
         break;
      n += syncsamp/4;
   }
   //fprintf(stderr, "Sync: %i\t%i\t%i\t%i\n", n, i, syncsamp/4, overlap);
   if ((syncsamp>>2)+i-SYNC_FIR < 0)
      die("Syncpoint missed by %i samples.", -((syncsamp>>2)+i-SYNC_FIR));
   // synced. From now no samples must get lost.
   fprintf(stderr, "Synced after %i samples. Read %i samples, Skip another %i samples\n", n+i, n+(syncsamp>>2), (syncsamp>>2)+i-SYNC_FIR);
   // discard until end of sync - overlap
   freadexact(inbuffer, 2*sizeof(short), (syncsamp>>2)+i-SYNC_FIR, in);

   // Scan !
   double fq = fmin;
   int findex = 0;
   while (fq < fmax)
   {  ++findex;
      int nfi = (int)(fq/freq*N + .5);
      if (findex < nfi)
         findex = nfi;
      // show
      fprintf(stderr, "Now at %.2f Hz ", findex * freq/N);
      // discard first samples
      freadexact(inbuffer, 2*sizeof(short), N, in);
      // write raw data #1
      FILE* fr = NULL;
      if (writeraw)
      {  char buf[30];
         sprintf(buf, "raw_%3f", findex * freq/N);
         buf[8] = 0;
         strcat(buf, ".dat");
         fr = fopen(buf, "w");
         write2ch(fr, inbuffer, N);
      }
      // show
      fputs("start...", stderr);
      // read data
      freadexact(inbuffer, 2*sizeof(short), N, in);
      // write raw data
      if (writeraw)
      {  write2ch(fr, inbuffer, N);
         fclose(fr);
         fr = NULL;
      }
      // show
      fputs("completed ", stderr);
      // analyze data
      analyze(findex);
      // show
      fprintf(stderr, "%f.1 dB\n", 20*log10(sqrt(sqr(ana[4]) + sqr(ana[5]))));
      // write result
      fprintf(resh, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%i\t%f\t%f\n", findex * freq/N,
       sqrt(sqr(ana[0]) + sqr(ana[1])), M_180_PI*atan2(ana[1], ana[0]), sqrt(sqr(ana[2]) + sqr(ana[3])), M_180_PI*atan2(ana[3], ana[2]), sqrt(sqr(ana[4]) + sqr(ana[5])), M_180_PI*atan2(ana[5], ana[4]),
       minmax[0], minmax[1], minmax[2], minmax[3], ana[6], ana[7]);
      fflush(resh);
      // next frequency
      fq *= fstep;
   }

}



// reference signal output
static void gensync()
{  short* dp = refbuffer;
   short* de = dp + syncsamp;
   while (dp != de)
   {  dp[0] = -32767;
      dp[1] = 32767;
      dp[2] = 0;
      dp[3] = 0;
      dp[4] = 32767;
      dp[5] = -32767;
      dp[6] = 0;
      dp[7] = 0;
      dp += 8;
   }
   // phase jump: Pi
   de = refbuffer + 2 * syncsamp;
   while (dp != de)
   {  dp[0] = 32767;
      dp[1] = -32767;
      dp[2] = 0;
      dp[3] = 0;
      dp[4] = -32767;
      dp[5] = 32767;
      dp[6] = 0;
      dp[7] = 0;
      dp += 8;
   }
}

static void genref(int fi)
{  // fill reference buffer
   double fs = M_2PI*fi/N;
   for (int i = 0; i < 2*N; i += 2)
      refbuffer[i+1] = -(refbuffer[i] = (short)(32767*cos(fs*i/2) + (double)rand()/(RAND_MAX+1)));
}

static void refplay(void*)
{  // we should increase the priority here
   #ifdef __OS2__
   DosSetPriority(PRTYS_THREAD, PRTYC_NOCHANGE, 1, 0);
   #else
   #error You need to implement a platform dependent way to increase the priority of this thread.
   #endif

   // write wav header
   _fsetmode(stdout, "b");
   static const char wavhdr[44] = {'R','I','F','F', -108,-1,-1,0x7f,
    'W','A','V','E','f','m','t',' ', 16,0,0,0, 1,0, 2,0, 0x80,0xbb,0,0, 0,0xee,2,0, 4,0, 16,0,
    'd','a','t','a', -144,-1,-1,0x7f};
   fwriteexact(wavhdr, 1, sizeof wavhdr, stdout);

   // pregap
   memset(refbuffer, 0, discardsamp * 2 * sizeof(short));
   fwriteexact(refbuffer, 2 * sizeof(short), discardsamp, stdout);
   // synchronize
   gensync();
   fwriteexact(refbuffer, 2 * sizeof(short), syncsamp, stdout);
   // overlap
   memset(refbuffer, 0, overlap * 2 * sizeof(short));
   fwriteexact(refbuffer, 2 * sizeof(short), overlap, stdout);
   // Wobble !
   double fq = fmin;
   int findex = 0;
   while (fq < fmax && !termrq)
   {  ++findex;
      int nfi = (int)(fq/freq*N + .5);
      if (findex < nfi)
         findex = nfi;
      genref(findex);
      // write reference
      fwriteexact(refbuffer, 2 * sizeof(short), N, stdout);
      // write reference 2nd try
      fwriteexact(refbuffer, 2 * sizeof(short), N, stdout);
      // next frequency
      fq *= fstep;
   }
}



/*static char* abbrev(const char* s, const char* token)
{  size_t l = strlen(s);
   return strnicmp(s, token, l) == 0 ? (char*)token + l : NULL;
}*/

static _cdecl void readint(const char* s, int* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%i%ln", r, &l) != 1 || l != strlen(s))
      die("Integer value expected, found %s", s);
}
static _cdecl void readintdef(const char* s, int* r, int d)
{  if (*s == 0)
      *r = d;
    else
   {  size_t l = (size_t)-1;
      if (sscanf(s, "%i%ln", r, &l) != 1 || l != strlen(s))
         die("Integer value expected, found %s", s);
}  }
static _cdecl void readfloat(const char* s, float* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%f%ln", r, &l) != 1 || l != strlen(s))
      die("Floating point value expected, found %s", s);
}
static _cdecl void readdouble(const char* s, double* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%lf%ln", r, &l) != 1 || l != strlen(s))
      die("Floating point value expected, found %s", s);
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
{  {"bn"  , (void(_cdecl*)(const char*,void*,int))&readN, &N, 0}
 , {"exec", (void(_cdecl*)(const char*,void*,int))&readstring, &execcmd, 0}
 , {"flog", (void(_cdecl*)(const char*,void*,int))&readdouble, &fstep, 0}
 , {"fmax", (void(_cdecl*)(const char*,void*,int))&readdouble, &fmax, 0}
 , {"fmin", (void(_cdecl*)(const char*,void*,int))&readdouble, &fmin, 0}
 , {"fq" ,  (void(_cdecl*)(const char*,void*,int))&readdouble, &freq, 0}
 , {"gd" ,  (void(_cdecl*)(const char*,void*,int))&setint, &gainmode, 3}
 , {"gg" ,  (void(_cdecl*)(const char*,void*,int))&setint, &gainmode, 2}
 , {"gr" ,  (void(_cdecl*)(const char*,void*,int))&setint, &gainmode, 1}
 , {"in" ,  (void(_cdecl*)(const char*,void*,int))&readstring, &infile, 0}
 , {"ln" ,  (void(_cdecl*)(const char*,void*,int))&readint, &loops, 1}
 , {"loop", (void(_cdecl*)(const char*,void*,int))&setint, &loops, INT_MAX}
 , {"ma" ,  (void(_cdecl*)(const char*,void*,int))&setint, &mode, 2}
 , {"mr" ,  (void(_cdecl*)(const char*,void*,int))&setint, &mode, 1}
 , {"psa",  (void(_cdecl*)(const char*,void*,int))&readint, &discardsamp, 0}
 , {"slvl", (void(_cdecl*)(const char*,void*,int))&readdouble, &synclevel, 0}
 , {"sov",  (void(_cdecl*)(const char*,void*,int))&readint, &overlap, 0}
 , {"sph",  (void(_cdecl*)(const char*,void*,int))&readdouble, &syncphase, 0}
 , {"sync", (void(_cdecl*)(const char*,void*,int))&readint, &syncsamp, 0}
 , {"v"  ,  (void(_cdecl*)(const char*,void*,int))&setflag, &verbose, true}
 , {"wd" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writedata, true}
 , {"wr" ,  (void(_cdecl*)(const char*,void*,int))&setflag, &writeraw, true}
 , {"zg" ,  (void(_cdecl*)(const char*,void*,int))&setint, &zeromode, 2}
 , {"zr" ,  (void(_cdecl*)(const char*,void*,int))&setint, &zeromode, 1}
};

static int searcharg(const char* arg, const char* elem)
{  return strnicmp(arg, elem, strlen(elem));
}

static void parsearg(const char* arg)
{  if (arg[0] == '@' || arg[0] == '<')
   {  // indirect file
      FILE* cf = fopen(arg+1, "r");
      if (cf == NULL)
         die("Failed to read command file %s.", arg+1);
      while (!feof(cf))
      {  char buffer[1024];
         fgets(buffer, sizeof buffer, cf);
         size_t l = strlen(buffer);
         if (l >= sizeof buffer-1 && buffer[sizeof buffer -2] != '\n')
            die("Line in command file exceeds 1023 characters.");
         // strip whitespaces
         while (strchr(" \t\r\n", buffer[--l]) != NULL)
            buffer[l] = 0;
         const char* ap = buffer;
         while (strchr(" \t\r\n", *ap) != NULL)
            ++ap;
         if (ap[0] == 0 || ap[0] == '#')
            continue; // skip empty lines and comments
         parsearg(ap); // THIS WILL NOT WORK WITH STRING ARGS !
      }
      return;
   }
   ArgMap* ap = (ArgMap*)bsearch(arg, argmap, sizeof argmap / sizeof *argmap, sizeof *argmap, (int (*)(const void*, const void*))&searcharg);
   if (ap == NULL)
      die("Illegal option %s.", arg);
   (*ap->func)(arg + strlen(ap->arg), ap->param, ap->iparam);
}

int main(int argc, char* argv[])
{  // parse cmdl
   while (--argc)
      parsearg(*++argv);

   syncsamp &= ~7; // must be multiple of 8

   /*resh = fdopen(3, "w");
   if (resh == NULL)
      die("Failed to open output handle 3 (%i).", errno);*/
   resh = stdout;

   switch (mode)
   {case 1:
      // generate
      refplay(NULL);
      break;
    case 3:
      // start reference generator
      _beginthread(refplay, NULL, 65536, NULL);
      sleep(1); // advantage sender, collect the data in the playrec buffer so far
    case 2:
      // analyze
      doanalysis(NULL);
      break;
    default:
      die("Unsupported mode: %i", mode);
   }

   // end reference player
   termrq = true;

   return 0;
}

