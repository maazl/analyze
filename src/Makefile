BINPATH = ../bin/

CC = gcc
CFLAGS = -std=c++11 -W -Wall -Wcast-qual -Wpointer-arith -Wcast-align -Wno-unused-result -fstrict-aliasing
#CFLAGS += -O6 -fomit-frame-pointer -malign-double
CFLAGS += -g
LDFLAGS = -lstdc++ -lm -lpthread
LDFLAGS += -s
LDFLAGSFFTW = -lfftw3
#EXEEXT = .exe


all : $(BINPATH)analyze$(EXEEXT) $(BINPATH)noise$(EXEEXT) $(BINPATH)ref$(EXEEXT) $(BINPATH)sweep$(EXEEXT) $(BINPATH)sampledist$(EXEEXT)

$(BINPATH)analyze$(EXEEXT) : analyze.cpp pca.h nrc.h nrc.cpp parser.cpp parser.h utils.cpp utils.h
	$(CC) $(CFLAGS) -o $@ analyze.cpp nrc.cpp parser.cpp utils.cpp $(LDFLAGSFFTW) $(LDFLAGS)

$(BINPATH)noise$(EXEEXT) : noise.cpp parser.cpp parser.h utils.cpp utils.h
	$(CC) $(CFLAGS) -o $@ noise.cpp parser.cpp utils.cpp $(LDFLAGSFFTW) $(LDFLAGS)

$(BINPATH)ref$(EXEEXT) : ref.cpp utils.cpp utils.h
	$(CC) $(CFLAGS) -o $@ ref.cpp utils.cpp $(LDFLAGSFFTW) $(LDFLAGS)

$(BINPATH)sampledist$(EXEEXT) : sampledist.cpp parser.cpp parser.h utils.cpp utils.h
	$(CC) $(CFLAGS) -o $@ sampledist.cpp parser.cpp utils.cpp $(LDFLAGS)

$(BINPATH)sweep$(EXEEXT) : sweep.cpp parser.cpp parser.h utils.cpp utils.h
	$(CC) $(CFLAGS) -o $@ sweep.cpp parser.cpp utils.cpp $(LDFLAGS)
