FFTWPATH = ../../k6fftwgel-1.1/

CFLAGS = -W -Wall -Zomf -Zcrtdll -Wcast-qual -Wpointer-arith -Wcast-align -pedantic
CFLAGS += -O6 -fomit-frame-pointer -malign-double -fstrict-aliasing -mcpu=k6
CFLAGS += -I$(FFTWPATH)fftw -I$(FFTWPATH)rfftw
#CFLAGS += -g
GCC3FIX = -lgcc3fix
LDFLAGS = -L$(FFTWPATH)fftw/.libs -L$(FFTWPATH)rfftw/.libs
LDFLAGS += -s -Zno-fork

BINPATH = ../bin/

ALL : $(BINPATH)Analyze.exe $(BINPATH)AnalyzeE.exe $(BINPATH)mls.exe $(BINPATH)ref.exe $(BINPATH)sweep.exe $(BINPATH)fft.exe $(BINPATH)sampledist.exe

$(BINPATH)Analyze.exe : Analyze.cpp pca.h nrc.h nrc.cpp
	gcc $(CFLAGS) -o $@ Analyze.cpp nrc.cpp $(LDFLAGS) -lrfftw -lfftw -lstdc++

$(BINPATH)AnalyzeE.exe : AnalyzeE.cpp
	gcc $(CFLAGS) -o $@ AnalyzeE.cpp $(LDFLAGS) -lrfftw -lfftw -lstdc++

$(BINPATH)mls.exe : mls.cpp
	gcc $(CFLAGS) -o $@ mls.cpp $(LDFLAGS) -lrfftw -lfftw -lstdc++

$(BINPATH)ref.exe : ref.cpp
	gcc $(CFLAGS) -o $@ ref.cpp $(LDFLAGS) -lrfftw -lfftw -lstdc++

$(BINPATH)fft.exe : fft.cpp
	gcc $(CFLAGS) -o $@ fft.cpp $(LDFLAGS) -lrfftw -lfftw -lstdc++

$(BINPATH)sampledist.exe : sampledist.cpp
	gcc $(CFLAGS) -o $@ sampledist.cpp $(LDFLAGS) -lstdc++

$(BINPATH)sweep.exe : sweep.cpp
	gcc $(CFLAGS) -o $@ sweep.cpp $(LDFLAGS)
