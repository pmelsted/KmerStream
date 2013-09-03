CC = g++
CXX = g++
INCLUDES = -I.
MAX_KMER_SIZE=64
CXXFLAGS = -c -Wall -Wno-reorder $(INCLUDES) -DMAX_KMER_SIZE=$(MAX_KMER_SIZE)
LDFLAGS =
LDLIBS  = -lm -lz

all: CXXFLAGS += -O3
all: target


debug: CXXFLAGS += -gstabs+ -O0
debug: LDFLAGS += -gstabs+
debug: target

profile: CXXFLAGS += -p -g -O2
profile: LDFLAGS += -p -g
profile: clean
profile: target

target: KmerStream


OBJECTS = lsb.o Kmer.o KmerIterator.o hash.o


KmerStream: KmerStream.o $(OBJECTS)
	$(CXX) $(INCLUDES) $(OBJECTS) KmerStream.o $(LDFLAGS) $(LDLIBS) -o KmerStream


lsb.o: lsb.cpp lsb.hpp
KmerStream.o: KmerStream.cpp StreamCounter.hpp
Kmer.o: Kmer.cpp Kmer.hpp
KmerIterator.o: KmerIterator.cpp KmerIterator.hpp
hash.o: hash.cpp hash.hpp

clean:
	rm -f *.o KmerStream
