FFTWPATH = ../../k6fftwgel-1.1/
BINPATH = ../bin/

CC = gcc
CFLAGS = -W -Wall -Zomf -Zcrtdll -Wcast-qual -Wpointer-arith -Wcast-align -pedantic
CFLAGS += -O6 -fomit-frame-pointer -malign-double -fstrict-aliasing -mcpu=k6
CFLAGS += -I$(FFTWPATH)fftw -I$(FFTWPATH)rfftw
#CFLAGS += -g
GCC3FIX = -lgcc3fix
LDFLAGS = -L$(FFTWPATH)fftw/.libs -L$(FFTWPATH)rfftw/.libs -lstdc++
LDFLAGS += -s -Zno-fork
LDFLAGSRFFTW = -lrfftw -lfftw
EXEEXT = .exe


ALL : $(BINPATH)Analyze$(EXEEXT) $(BINPATH)AnalyzeE$(EXEEXT) $(BINPATH)mls$(EXEEXT) $(BINPATH)ref$(EXEEXT) $(BINPATH)sweep$(EXEEXT) $(BINPATH)fft$(EXEEXT) $(BINPATH)sampledist$(EXEEXT)

$(BINPATH)Analyze$(EXEEXT) : Analyze.cpp pca.h nrc.h nrc.cpp
	$(CC) $(CFLAGS) -o $@ Analyze.cpp nrc.cpp $(LDFLAGSRFFTW) $(LDFLAGS)

$(BINPATH)AnalyzeE$(EXEEXT) : AnalyzeE.cpp
	$(CC) $(CFLAGS) -o $@ AnalyzeE.cpp $(LDFLAGSRFFTW) $(LDFLAGS)

$(BINPATH)noise$(EXEEXT) : noise.cpp
	$(CC) $(CFLAGS) -o $@ mls.cpp $(LDFLAGSRFFTW) $(LDFLAGS)

$(BINPATH)ref$(EXEEXT) : ref.cpp
	$(CC) $(CFLAGS) -o $@ ref.cpp $(LDFLAGSRFFTW) $(LDFLAGS)

$(BINPATH)fft$(EXEEXT) : fft.cpp
	$(CC) $(CFLAGS) -o $@ fft.cpp $(LDFLAGSRFFTW) $(LDFLAGS)

$(BINPATH)sampledist$(EXEEXT) : sampledist.cpp
	$(CC) $(CFLAGS) -o $@ sampledist.cpp $(LDFLAGS)

$(BINPATH)sweep$(EXEEXT) : sweep.cpp
	$(CC) $(CFLAGS) -o $@ sweep.cpp $(LDFLAGS)
