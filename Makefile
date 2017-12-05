CC = gcc
CXX = clang++
CFLAGS = -Wall -Wextra -Wshadow -Wformat-nonliteral -Wformat-security -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -O2 -O3 
CXXFLAGS = -std=c++11
LDLIBS = -lm -larmadillo -lblas -llapack -lfftw3 libz.a -lboost_system -lboost_filesystem

programs = edf2cfs

all: $(programs)

$(programs): SHA1.o edflib.o resample.o

clean:
	$(RM) *.o $(programs) *.[be]df
