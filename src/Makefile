FFTWPATH = ../../fftw-2.1.5/
BINPATH = ../bin/

CC = gcc
CFLAGS = -W -Wall -Wcast-qual -Wpointer-arith -Wcast-align -pedantic
CFLAGS += -O6 -fomit-frame-pointer -malign-double -fstrict-aliasing
CFLAGS += -I$(FFTWPATH)fftw -I$(FFTWPATH)rfftw
#CFLAGS += -g
#GCC3FIX = -lgcc3fix
LDFLAGS = -L$(FFTWPATH)fftw/.libs -L$(FFTWPATH)rfftw/.libs -lstdc++
LDFLAGS += -s
LDFLAGSRFFTW = -lrfftw -lfftw
#EXEEXT = .exe


ALL : $(BINPATH)analyze$(EXEEXT) $(BINPATH)noise$(EXEEXT) $(BINPATH)ref$(EXEEXT) $(BINPATH)sweep$(EXEEXT) $(BINPATH)sampledist$(EXEEXT)

$(BINPATH)analyze$(EXEEXT) : analyze.cpp pca.h nrc.h nrc.cpp parser.cpp parser.h
	$(CC) $(CFLAGS) -o $@ analyze.cpp nrc.cpp parser.cpp $(LDFLAGSRFFTW) $(LDFLAGS)

$(BINPATH)analyzeE$(EXEEXT) : analyzeE.cpp
	$(CC) $(CFLAGS) -o $@ analyzeE.cpp $(LDFLAGSRFFTW) $(LDFLAGS)

$(BINPATH)noise$(EXEEXT) : noise.cpp parser.cpp parser.h
	$(CC) $(CFLAGS) -o $@ noise.cpp parser.cpp $(LDFLAGSRFFTW) $(LDFLAGS)

$(BINPATH)ref$(EXEEXT) : ref.cpp
	$(CC) $(CFLAGS) -o $@ ref.cpp $(LDFLAGSRFFTW) $(LDFLAGS)

$(BINPATH)sampledist$(EXEEXT) : sampledist.cpp
	$(CC) $(CFLAGS) -o $@ sampledist.cpp $(LDFLAGS)

$(BINPATH)sweep$(EXEEXT) : sweep.cpp parser.cpp parser.h
	$(CC) $(CFLAGS) -o $@ sweep.cpp parser.cpp $(LDFLAGS)
