#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#include <rfftw.h>


static const int fsamp = 48000; // well, not nice...

void die(int rc, const char* msg, ...)
{  va_list va;
   va_start(va, msg);
   vfprintf(stderr, msg, va);
   va_end(va);
   exit(rc);
}

void wavwriter(FILE* fo, size_t nsamp)
{  // fake wav header
   int wavhdr[11] =
   /*44100
   { 0x46464952, 0xFF, 0xFF, 0xFF, 0x7F, 0x57, 0x41, 0x56, 0x45, 0x66, 0x6D, 0x74, 0x20,
     0x10, 0x00, 0x00, 0x00, 0x01, 0x00, 0x01, 0x00, 0x44, 0xAC, 0x00, 0x00, 0x88, 0x58, 0x01, 0x00,
     0x02, 0x00, 0x10, 0x00, 0x64, 0x61, 0x74, 0x61, 0xFF, 0xFF, 0xFF, 0x7F };*/
   // 48000
   { 0x46464952, 0xffffffff, 0x45564157, 0x20746D66,
     0x00000010, 0x00020001, 0x0000BB80, 0x0002ee00,
     0x00100004, 0x61746164, 0xffffffff };

   wavhdr[10] = nsamp * sizeof(short);
   wavhdr[1] = wavhdr[10] + 44;

   fwrite(wavhdr, sizeof wavhdr, 1, fo);
}

inline static double myrand()
{  unsigned long l = (((rand() << 11) ^ rand()) << 11) ^ rand();
   return (double)l / ULONG_MAX;
}


static double equalscale(int)
{  return 1;
}

static double f_scale(int i)
{  return i;
}

static double sqrt_f_scale(int i)
{  return sqrt(i);
}

static double scalepow = 0;

static double powerscale(int i)
{  return pow(i, scalepow);
}

static double randomphase(int)
{  return 2*M_PI * myrand();
}

static double phaseparam;

static double constphase(int)
{  return phaseparam;
}

static double (*scalefn)(int i) = equalscale;
static double (*phasefn)(int i) = randomphase;


void _cdecl setscalefunc(const char* name, double (**scalefn)(int i))
{  static const struct scfne
   {  char     name[8];
      double   (*fn)(int);
   } scfnt[] =
   {  {"eq",     &equalscale}
    , {"f",      &f_scale}
    , {"sqrt_f", &sqrt_f_scale}
   };

   if (name == NULL || *name == 0)
      return;
   const scfne* sp = (const scfne*)bsearch(name, scfnt, sizeof(scfnt)/sizeof(*scfnt), sizeof(*scfnt), (int (*)(const void*, const void*))&stricmp);
   if (sp == NULL)
   {  //die (34, "invalid scale function %s", name);
      *scalefn = &powerscale;
      scalepow = atof(name);
      return;
   }
   *scalefn = sp->fn;
}

void _cdecl setphasefunc(const char* name, double (**phasefn)(int i))
{  static const struct phfne
   {  char     name[8];
      double   (*fn)(int);
   } phfnt[] =
   {  {"rand",   &randomphase}
   };

   if (name == NULL || *name == 0)
      return;
   const phfne* sp = (const phfne*)bsearch(name, phfnt, sizeof(phfnt)/sizeof(*phfnt), sizeof(*phfnt), (int (*)(const void*, const void*))&stricmp);
   if (sp == NULL)
   {  //die (34, "invalid scale function %s", name);
      *phasefn = &constphase;
      phaseparam = atof(name);
      return;
   }
   *phasefn = sp->fn;
}

int main(int argc, char**argv)
{  srand(clock());

   // command line
   const char* dfile = NULL;
   size_t nrep = 1;
   switch (argc)
   {default:
      die(45, "usage: ref nsamp fmin fmax [scalefn[,phasefn] [nrep [file]]]");
    case 7:
      dfile = argv[6];
    case 6:
      nrep = atol(argv[5]);
    case 5:
      setscalefunc(strtok(argv[4], "_,;"), &scalefn);
      setphasefunc(strtok(NULL, ""), &phasefn);
    case 4:;
   }
   size_t nsamp = atol(argv[1]);
   double fmin = atof(argv[2]);
   double fmax = atof(argv[3]);

   // round fmin & fmax
   double finc = (double)fsamp / nsamp;
   fprintf(stderr, "finc = %g\n", finc);

   int imin = (int)floor((fmin / finc) +.5);
   int imax = (int)floor((fmax / finc) +.5);
   fprintf(stderr, "fmin = %g\n"
                   "fmax = %g\n", imin*finc, imax*finc);
   if ((unsigned int)imax > nsamp/2 || (unsigned int)imin > (unsigned int)imax)
      die(34, "fmin and/or fmax out of range");

   // generate coefficients
   float* fftbuf = new float[nsamp+1];
   if (fftbuf == NULL)
      die(39, "malloc(%lu) failed", nsamp);
   memset(fftbuf, 0, sizeof fftbuf);   // all coefficients -> 0
   for (int i = imin; i <= imax; ++i)
   {  double r = (*scalefn)(i);
      double phi = (*phasefn)(i);
      fftbuf[nsamp-i] = r * sin(phi);  // b[i]
      fftbuf[i] = r * cos(phi);        // a[i]
   }

   // ifft
   float* sampbuf = new float[nsamp];
   rfftw_plan plan = rfftw_create_plan(nsamp, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
   rfftw_one(plan, fftbuf, sampbuf);   // IFFT
   rfftw_destroy_plan(plan);
   delete[] fftbuf; // no longer needed

   // normalize
   double fnorm = 0;
   const float* sp = sampbuf;
   const float* spe = sp + nsamp;
   for (; sp != spe; ++sp)
      if (fabs(*sp) > fnorm)
         fnorm = fabs(*sp);
   fnorm = 32767/fnorm;
   // and quantize
   short* buf = new short[2*nsamp];
   sp = sampbuf;
   short* dp = buf;
   while (sp != spe)
   {  register short s = (short)floor(*sp++ * fnorm + myrand());
      dp[0] = s;
      dp[1] = -s;
      dp += 2;
   }
   delete[] sampbuf; // no longer needed

   // write
   FILE* of;
   if (dfile)
   {  of = fopen(dfile, "wb");
      if (of == NULL)
         die (41, "Failed to open %s for writing", dfile);
      wavwriter(of, 2*nsamp * nrep);
   } else // stdout
   {  _fsetmode(stdout, "b");
      of = stdout;
      wavwriter(of, 0x1fffffdc);
   }
   while (--nrep)
      fwrite(buf, 2*nsamp * sizeof(short), 1, of);

   // cleanup, well not in case of an exception...
   delete[] buf;
}

