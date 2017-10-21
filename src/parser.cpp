#include "parser.h"

#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <process.h>

bool termrq = false;

void die(int rc, const char* msg, ...)
{  termrq = true;
   va_list va;
   va_start(va, msg);
   vfprintf(stderr, msg, va);
   va_end(va);
   fputc('\n', stderr);
   exit(rc);
}

static int searcharg(const char* arg, const char* elem)
{  return strnicmp(arg, elem, strlen(elem));
}

void parsearg(const char* arg)
{  if (arg[0] == '@' || arg[0] == '<')
   {  // indirect file
      FILE* cf = fopen(arg+1, "r");
      if (cf == NULL)
         die(37, "Failed to read command file %s.", arg+1);
      while (!feof(cf))
      {  char buffer[1024];
         fgets(buffer, sizeof buffer, cf);
         size_t l = strlen(buffer);
         if (l >= sizeof buffer-1 && buffer[sizeof buffer -2] != '\n')
            die(39, "Line in command file exceeds 1023 characters.");
         // strip whitespaces
         while (strchr(" \t\r\n", buffer[--l]) != NULL)
            buffer[l] = 0;
         if (buffer[0] == 0)
            continue; // skip empty lines
         const char* ap = buffer;
         while (strchr(" \t\r\n", *ap) != NULL)
            ++ap;
         if (ap[0] == '#')
            continue; // skip comments
         parsearg(ap); // THIS WILL NOT WORK WITH STRING ARGS !
      }
      return;
   }
   ArgMap* ap = (ArgMap*)bsearch(arg, argmap, argmap_size, sizeof *argmap, (int (*)(const void*, const void*))&searcharg);
   if (ap == NULL)
      die(44, "Illegal option %s.", arg);
   arg += strlen(ap->arg);
   if (*arg == '=')
      ++arg;
   (*ap->func)(arg, ap->param, ap->iparam);
}

void readint(const char* s, int* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%i%zn", r, &l) != 1 || l != strlen(s))
      die(42, "Integer value expected, found %s", s);
}
void readintdef(const char* s, int* r, int d)
{  if (*s == 0)
      *r = d;
    else
   {  size_t l = (size_t)-1;
      if (sscanf(s, "%i%zn", r, &l) != 1 || l != strlen(s))
         die(42, "Integer value expected, found %s", s);
}  }
void readuint(const char* s, unsigned int* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%u%zn", r, &l) != 1 || l != strlen(s))
      die(42, "Unsigned integer value expected, found %s", s);
}
void readuintdef(const char* s, unsigned int* r, unsigned int d)
{  if (*s == 0)
      *r = d;
    else
   {  size_t l = (size_t)-1;
      if (sscanf(s, "%u%n", r, &l) != 1 || l != strlen(s))
         die(42, "Unsigned integer value expected, found %s", s);
}  }
void readfloat(const char* s, float* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%f%zn", r, &l) != 1 || l != strlen(s))
      die(42, "Floating point value expected, found %s", s);
}
void readdouble(const char* s, double* r)
{  size_t l = (size_t)-1;
   if (sscanf(s, "%lf%zn", r, &l) != 1 || l != strlen(s))
      die(42, "Floating point value expected, found %s", s);
}

void readstring(const char* s, const char** cpp)
{  *cpp = strdup(s);
}
void readstringdef(const char* s, const char** cpp, const char* def)
{  *cpp = strdup(*s ? s : def);
}

void setflag(const char* s, bool* r)
{  if (*s == 0 || strcmp(s, "1") == 0 || stricmp(s, "true") == 0)
      *r = true;
   else if (strcmp(s, "0") == 0 || stricmp(s, "false") == 0)
      *r = false;
   else
      die(42, "Invalid boolean value %s", s);
}

void setint(const char* s, int* r, int v)
{  if (*s)
      die(42, "Option does not have parameters");
   *r = v;
}

void setbit(const char* s, unsigned int* r, unsigned int v)
{  if (*s)
      die(42, "Option does not have parameters");
   *r |= v;
}


