BINPATH = ../bin/

CC = gcc
CFLAGS = -std=c++11 -W -Wall -Wcast-qual -Wpointer-arith -Wcast-align -Wno-unused-result
CFLAGS += -O6 -fomit-frame-pointer -malign-double -fstrict-aliasing
#CFLAGS += -g
#GCC3FIX = -lgcc3fix
LDFLAGS = -lstdc++ -lm -lpthread
LDFLAGS += -s
LDFLAGSFFTW = -lfftw3
#EXEEXT = .exe


all : $(BINPATH)analyze$(EXEEXT) $(BINPATH)noise$(EXEEXT) $(BINPATH)ref$(EXEEXT) $(BINPATH)sweep$(EXEEXT) $(BINPATH)sampledist$(EXEEXT)

$(BINPATH)analyze$(EXEEXT) : analyze.cpp pca.h nrc.h nrc.cpp parser.cpp parser.h
	$(CC) $(CFLAGS) -o $@ analyze.cpp nrc.cpp parser.cpp $(LDFLAGSFFTW) $(LDFLAGS)

$(BINPATH)noise$(EXEEXT) : noise.cpp parser.cpp parser.h
	$(CC) $(CFLAGS) -o $@ noise.cpp parser.cpp $(LDFLAGSFFTW) $(LDFLAGS)

$(BINPATH)ref$(EXEEXT) : ref.cpp
	$(CC) $(CFLAGS) -o $@ ref.cpp $(LDFLAGSFFTW) $(LDFLAGS)

$(BINPATH)sampledist$(EXEEXT) : sampledist.cpp parser.cpp parser.h
	$(CC) $(CFLAGS) -o $@ sampledist.cpp parser.cpp $(LDFLAGS)

$(BINPATH)sweep$(EXEEXT) : sweep.cpp parser.cpp parser.h
	$(CC) $(CFLAGS) -o $@ sweep.cpp parser.cpp $(LDFLAGS)
