#include <stdlib.h>

// Termination flag, set by die(...).
extern bool termrq;

// Terminate the application with an error message and set termrq
// to notify other threads (if any).
void die(int rc, const char* msg, ...);

typedef void (*ArgFn)(const char* rem, void* param, int iparam);

// Dispatch table for command line argument lookup
extern const struct ArgMap
{  char  arg[8];
   ArgFn func;
   void* param;
   int   iparam;
} argmap[];
// Number of entries in argmap[].
// Must be set to sizeof argmap / sizeof *argmap once argmap is a complete type.
extern const size_t argmap_size;

// Parse the specified argument.
// The argument may start with a prefix from the dispatch table
// or have the format @filename in which case parsearg is called for each line
// in the file recursively.
void parsearg(const char* arg);

// Helper functions for the Dispatch table ...

// Read integer parameter
// *param = &intparam
void readint(const char* s, int* r);
// Read integer parameter with default value
// *param = &intparam
// iparam = default value if argument is specified without a value.
void readintdef(const char* s, int* r, int d);
// Read integer parameter
// *param = &intparam
void readuint(const char* s, unsigned int* r);
// Read integer parameter with default value
// *param = &intparam
// iparam = default value if argument is specified without a value.
void readuintdef(const char* s, unsigned int* r, unsigned int d);
// Read float parameter 
// *param = &floatparam
void readfloat(const char* s, float* r);
// Read double parameter 
// *param = &doubleparam
void readdouble(const char* s, double* r);

// Read string parameter
// *param = &charpointer
// This function returns a copy of the argument with strdup that should
// be freed with free().
void readstring(const char* s, const char** cpp);
// Read string parameter with default value
// *param = &charpointer
// iparam = "default value"
// This function returns a copy of the argument with strdup that should
// be freed with free().
void readstringdef(const char* s, const char** cpp, const char* def);

// Set boolean parameter
// *param = &bool
// This has no arguments. It unconditional sets the parameter to true.
void setflag(const char* s, bool* r);

// Set integer parameter to a constant value
// *param = &intparam
// iparam = constant
// This has no arguments. The constant is unconditionally assigned to *param.
// Be careful with enum types because their size may differ. 
void setint(const char* s, int* r, int v);
#define setuint setint

// Set a costant bit in an integer parameter
// *param = &intparam
// ivalue = bit_number
// This has no arguments. The n-th bit in *param ios set.
void setbit(const char* s, unsigned int* r, unsigned int v);

