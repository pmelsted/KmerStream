INCLUDES = -I. 
MAX_KMER_SIZE=64
CXXFLAGS = -c -pthread -Wall -Wno-reorder -Wno-header-guard $(INCLUDES) -DMAX_KMER_SIZE=$(MAX_KMER_SIZE) 
LDFLAGS =


OPENMP=
#ifeq ($(findstring clang,$(shell $(CXX) --version)), clang) 
#OPENMP=
#else
#OPENMP=-lgomp
#CXXFLAGS += -fopenmp
#endif

LDLIBS  = -lm -lz $(OPENMP) -lpthread


all: CXXFLAGS += -O3
all: target



debug: CXXFLAGS += -g -O0
debug: LDFLAGS += -g
debug: target

profile: CXXFLAGS += -p -g -O2
profile: LDFLAGS += -p -g
profile: clean
profile: target

target: KmerStream


OBJECTS = lsb.o Kmer.o KmerIterator.o hash.o RepHash.o


KmerStream: KmerStream.o $(OBJECTS)
	$(CXX) $(INCLUDES) $(OBJECTS) KmerStream.o $(LDFLAGS) $(LDLIBS) -o KmerStream


lsb.o: lsb.cpp lsb.hpp
KmerStream.o: KmerStream.cpp StreamCounter.hpp
Kmer.o: Kmer.cpp Kmer.hpp
KmerIterator.o: KmerIterator.cpp KmerIterator.hpp
hash.o: hash.cpp hash.hpp
RepHash.o: RepHash.cpp RepHash.hpp

clean:
	rm -f *.o KmerStream
