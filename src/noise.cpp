#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#include <rfftw.h>

#include <complex>
using namespace std;
typedef complex<double> Complex;

#include "parser.h"

#define M_180_PI (180./M_PI)


static unsigned n_fft = 0;
static unsigned f_samp = 48000; // sampling rate
static double f_min = 0;
static double f_max = 0;
static double f_inc = 1;
static double f_log = 0;
static double scalepow = 0;
static unsigned n_harmonic = 0;
static bool stereo = false;
static const char* F_data = NULL;
static const char* F_res = NULL;
static const char* F_wav = NULL;
static unsigned n_rep = 1;
//static int chirpphase = 0;

void wavwriter(FILE* fo, size_t n_samp)
{  // fake wav header
   int wavhdr[11] =
   // 48000
   { 0x46464952, 0xffffffff, 0x45564157, 0x20746D66,
     0x00000010, 0x00020001, 0x0000BB80, 0x0002ee00,
     0x00100004, 0x61746164, 0xffffffff };

   // patch data size
   wavhdr[10] = n_samp * sizeof(short);
   wavhdr[1] = wavhdr[10] + 44;
   // patch sampling rate
   wavhdr[6] = f_samp;
   wavhdr[7] = f_samp << 2; // * 4 bytes per sample

   fwrite(wavhdr, sizeof wavhdr, 1, fo);
}

inline static double myrand()
{  unsigned long l = (((rand() << 11) ^ rand()) << 11) ^ rand();
   return (double)l / ULONG_MAX;
}


static void readN(const char* s, unsigned* r)
{  bool ex;
   if (ex = *s == '^')
      ++s;
   readuint(s, r);
   if (ex)
      *r = 1 << *r;
}

const ArgMap argmap[] = // must be sorted
{  {"bn"  , (ArgFn)&readN, &n_fft, 0}
 //, {"chirp",(ArgFn)&readintdef, &chirpphase, 1}
 , {"finc", (ArgFn)&readdouble, &f_inc, 0}
 , {"flog", (ArgFn)&readdouble, &f_log, 0}
 , {"fmax", (ArgFn)&readdouble, &f_max, 0}
 , {"fmin", (ArgFn)&readdouble, &f_min, 0}
 , {"fsamp",(ArgFn)&readuint,   &f_samp, 0}
 , {"harm", (ArgFn)&readuintdef,&n_harmonic, 3}
 , {"ln" ,  (ArgFn)&readuint,   &n_rep, 1}
 , {"loop", (ArgFn)&setuint,    &n_rep, INT_MAX}
 , {"mst",  (ArgFn)&setflag,    &stereo, true}
 , {"scale",(ArgFn)&readdouble, &scalepow, 0}
 , {"wd",   (ArgFn)&readstring, &F_data, 0}
 , {"wr",   (ArgFn)&readstring, &F_res, 0}
 , {"ww",   (ArgFn)&readstring, &F_wav, 0}
};
const size_t argmap_size = sizeof argmap / sizeof *argmap;


int main(int argc, char**argv)
{  srand(clock());

   // parse cmdl
   while (--argc)
      parsearg(*++argv);

   if (n_fft == 0)
      die(48, "usage: noise bn<fft_size> fmin<min_freq> fmax<max_freq> [wd<data_file>] [wr<wave_file>]\n"
              "See documentation for more options.");

   // round fmin & fmax
   double f_bin = (double)f_samp/n_fft;
   size_t i_min = (int)floor(f_min/f_bin +.5);
   size_t i_max = (int)floor(f_max/f_bin +.5);
   if (i_max > n_fft/2 || i_min > i_max)
      die(34, "fmin and/or fmax out of range");
   f_inc -= .5; // increment compensation + rounding
   f_log += 1;

   // generate coefficients
   float* fftbuf = new float[(stereo+1)*(n_fft+1)];
   int* harmonics = new int[n_fft+2];
   if (fftbuf == NULL || harmonics == NULL)
      die(39, "malloc(%lu) failed", n_fft);
   memset(fftbuf, 0, (stereo+1)*(n_fft+1) * sizeof *fftbuf);   // all coefficients -> 0
   memset(harmonics, 0, (n_fft/2+1) * sizeof *harmonics);
   int sign = 1;
   for (size_t i = i_min; i <= i_max; ++i)
   {  for (size_t j = 1; j <= n_harmonic && i*j <= n_fft/2; ++j)
         if (harmonics[i*j])
            goto next_f;
      for (size_t j = 1; i*j <= n_fft/2; ++j)
         harmonics[i*j] = j * sign;
      {  double r = pow(i, scalepow);
         double phi = 2*M_PI * myrand();
         fftbuf[n_fft-i] = r * sin(phi);  // b[i]
         fftbuf[i] = r * cos(phi);        // a[i]
         if (sign < 0)
         {  fftbuf[2*n_fft+1-i] = fftbuf[n_fft-i];
            fftbuf[n_fft+1+i] = fftbuf[i];
         }
      }
      if (stereo)
         sign = -sign;
      i = (int)floor(i * f_log + f_inc);
    next_f:;
   }

   // write design result
   if (F_data)
   {  FILE* of = fopen(F_data, "w");
      if (of == NULL)
         die (41, "Failed to open %s for writing", F_data);
      for (size_t i = 0; i <= n_fft/2; ++i)
      {  Complex ci(fftbuf[i], i && i != n_fft/2 ? fftbuf[n_fft-i] : 0);
         //       freq |A|  arg A ai  bi
         fprintf(of, "%12g %12g %12g %12g %12g %8i\n",
            i*f_bin, abs(ci), arg(ci)*M_180_PI, ci.real(), ci.imag(), harmonics[i]); 
      }   
      fclose(of);
   }

   if (!F_wav && !F_res)
      goto end;

   if (!stereo)
   {
      // ifft
      float* sampbuf = new float[n_fft];
      rfftw_plan plan = rfftw_create_plan(n_fft, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
      rfftw_one(plan, fftbuf, sampbuf);   // IFFT
      rfftw_destroy_plan(plan);

      // normalize
      double fnorm = 0;
      float* sp = sampbuf;
      const float* const spe = sp + n_fft;
      for (; sp != spe; ++sp)
         if (fabs(*sp) > fnorm)
            fnorm = fabs(*sp);
      fnorm = 1/fnorm;
      for (sp = sampbuf; sp != spe; ++sp)
         *sp *= fnorm; 
         
      // write result
      if (F_res)
      {  FILE* of = fopen(F_res, "w");
         if (of == NULL)
            die (41, "Failed to open %s for writing", F_res);
         for (sp = sampbuf; sp != spe; ++sp)
            fprintf(of, "%12g\n", *sp);
         fclose(of);
      }   
         
      // and quantize
      if (F_wav)
      {  short* buf = new short[2*n_fft];
         sp = sampbuf;
         short* dp = buf;
         while (sp != spe)
         {  register short s = (short)floor(*sp++ * 32767 + myrand());
            dp[0] = s;
            dp[1] = -s;
            dp += 2;
         }
         delete[] sampbuf; // no longer needed

         FILE* of;
         if (strcmp(F_wav, "-") != 0)
         {  of = fopen(F_wav, "wb");
            if (of == NULL)
               die (41, "Failed to open %s for writing", F_wav);
            wavwriter(of, 2*n_fft * n_rep);
         } else // stdout
         {  _fsetmode(stdout, "b");
            of = stdout;
            wavwriter(of, 0x1fffffdc);
         }
         do fwrite(buf, 2*n_fft * sizeof(short), 1, of);
         while (--n_rep);

         // cleanup, well not in case of an exception...
         delete[] buf;
      }

   } else // stereo
   {
      // split channels
      float* const fftr = fftbuf + n_fft+1;
      for (size_t i = i_min; i <= i_max; ++i)
      {  if (harmonics[i] < 0)
            fftbuf[i] = fftbuf[n_fft-i] = 0;
      }
   
      // ifft
      float* sampbuf = new float[2*n_fft];
      rfftw_plan plan = rfftw_create_plan(n_fft, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
      rfftw_one(plan, fftbuf, sampbuf);   // IFFT
      rfftw_one(plan, fftr, sampbuf+n_fft);
      rfftw_destroy_plan(plan);

      // normalize
      {  double fnorm = 0;
         float* sp = sampbuf;
         const float* const spe = sp + 2*n_fft;
         for (; sp != spe; ++sp)
            if (fabs(*sp) > fnorm)
               fnorm = fabs(*sp);
         fnorm = 1/fnorm;
         for (sp = sampbuf; sp != spe; ++sp)
            *sp *= fnorm;
      } 
         
      // write result
      if (F_res)
      {  FILE* of = fopen(F_res, "w");
         if (of == NULL)
            die (41, "Failed to open %s for writing", F_res);
         const float* const spe = sampbuf + n_fft;
         for (float* sp = sampbuf; sp != spe; ++sp)
            fprintf(of, "%12g %12g %12g\n", sp[0]+sp[n_fft], sp[0], sp[n_fft]);
         fclose(of);
      }   
         
      // and quantize
      if (F_wav)
      {  short* buf = new short[2*n_fft];
         float* sp = sampbuf;
         const float* const spe = sp + n_fft;
         short* dp = buf;
         while (sp != spe)
         {  dp[0] = (short)floor(sp[0] * 32767 + myrand());
            dp[1] = (short)floor(sp[n_fft] * 32767 + myrand());
            ++sp;
            dp += 2;
         }
         delete[] sampbuf; // no longer needed

         FILE* of;
         if (strcmp(F_wav, "-") != 0)
         {  of = fopen(F_wav, "wb");
            if (of == NULL)
               die (41, "Failed to open %s for writing", F_wav);
            wavwriter(of, 2*n_fft * n_rep);
         } else // stdout
         {  _fsetmode(stdout, "b");
            of = stdout;
            wavwriter(of, 0x1fffffdc);
         }
         do fwrite(buf, 2*n_fft * sizeof(short), 1, of);
         while (--n_rep);

         // cleanup, well not in case of an exception...
         delete[] buf;
      }
   }
   
 end:
   delete[] harmonics;
   delete[] fftbuf;
}

