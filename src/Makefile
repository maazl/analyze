CFLAGS = -W -Wall -Zomf -Zcrtdll -Wcast-qual -Wpointer-arith -Wcast-align -pedantic
CFLAGS += -O6 -fomit-frame-pointer -malign-double -fstrict-aliasing -mcpu=k6
CFLAGS += -I../k6fftwgel-1.1/fftw -I../k6fftwgel-1.1/rfftw
#CFLAGS += -g
GCC3FIX = -lgcc3fix
LDFLAGS = -L../k6fftwgel-1.1/fftw/.libs -L../k6fftwgel-1.1/rfftw/.libs
LDFLAGS += -s

BINPATH = ../bin/

ALL : $(BINPATH)Analyze.exe $(BINPATH)ref.exe $(BINPATH)sweep.exe $(BINPATH)fft.exe $(BINPATH)sampledist.exe

$(BINPATH)Analyze.exe : Analyze.cpp pca.h nrc.h nrc.cpp makefile
	gcc $(CFLAGS) -o $@ Analyze.cpp nrc.cpp $(LDFLAGS) -lrfftw -lfftw -lstdc++

$(BINPATH)AnalyzeE.exe : AnalyzeE.cpp makefile
	gcc $(CFLAGS) -o $@ AnalyzeE.cpp $(LDFLAGS) -lrfftw -lfftw -lstdc++

$(BINPATH)ref.exe : ref.cpp makefile
	gcc $(CFLAGS) -o $@ ref.cpp $(LDFLAGS) -lrfftw -lfftw -lstdc++

$(BINPATH)fft.exe : fft.cpp makefile
	gcc $(CFLAGS) -o $@ fft.cpp $(LDFLAGS) -lrfftw -lfftw -lstdc++

$(BINPATH)sampledist.exe : sampledist.cpp makefile
	gcc $(CFLAGS) -o $@ sampledist.cpp $(LDFLAGS) -lstdc++

$(BINPATH)sweep.exe : sweep.cpp
	gcc $(CFLAGS) -o $@ sweep.cpp $(LDFLAGS)

$(BINPATH)AnalyzeE.exe : AnalyzeE.cpp makefile
	gcc $(CFLAGS) -o $@ AnalyzeE.cpp $(LDFLAGS) -lrfftw -lfftw
