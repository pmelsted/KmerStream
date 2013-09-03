#ifndef STREAMCOUNTER_HPP
#define STREAMCOUNTER_HPP

#include <sstream>
#include <stdint.h>
#include <string>
#include <cmath>

#include "lsb.hpp"

class StreamCounter{
public:
  StreamCounter(double e_, size_t seed_) : MAX_TABLE(32), maxVal(3ULL), e(e_), seed(seed_) {
    size_t numbits = (size_t)(48.0/(e*e) + 1); // approx 3 std-dev, true with 0.001 prob.
    if (numbits < 8192) { numbits = 8192;}

    size =(numbits+31)/32; 
    size--; size |= size>>1; size |= size>>2; size |= size>>4; size |= size>>8; size |= size>>16; size |= size>>32; size++;
    mask = (size<<5) -1; // 2^5=32  per block
    M = new size_t[MAX_TABLE];
    memset(M,0,MAX_TABLE);
    table = new uint64_t[size*MAX_TABLE];
    memset(table, 0, size*MAX_TABLE);
  }

  ~StreamCounter() { 
    delete[] table;
    delete[] M;
  }
  
    
  void operator()(uint64_t hashval) {
    // hashval is XXX .. XXX1000.. (w-1 times) ..00 0
    size_t w = bitScanForward(hashval); // 1-based index of first 1
    if (w >= MAX_TABLE) {
      w = MAX_TABLE-1;
    }
    if (M[w] == size*32*maxVal) { // ok full
      return;
    }
    // hashval is now XXXX...XX, random
    hashval = hashval >> (w+1); // shift away pattern of XXX10000
    
    uint64_t index = hashval & mask;
    uint64_t val = getVal(index,w);
    if (val != maxVal) { // max count
      setVal(index,w,val+1);
      M[w]++;
    }


  }

  bool join(const StreamCounter &o) {
    if (size != o.size || seed != o.seed) {
      return false;
    }

    for (size_t i = 0; i < MAX_TABLE; i++) {
      for (size_t j = 0; j < size*32; j++) {
	setVal(j,i,getVal(j,i)+o.getVal(j,i));
      }
    }

    return true;
  }
  
  size_t F0() const {
    size_t R = size*32;
    double sum = 0;

    int n = 0;
    double limit = 0.2;
    while (n == 0) {
      for (size_t i = 0; i < MAX_TABLE; i++) {
	size_t ts = 0;
	for (size_t j = 0; j < R; j++) {
	  if (getVal(j,i) > 0) {
	    ts++;
	  }
	}
	
	if (ts <= (1-limit)*R && ts >= limit*R) {
	  double est = (log(1.0-ts/((double) R))/log(1.0-1.0/R)) * pow(2.0,i);
	  std::cout << i << " " << ts << " " << limit<<  " " << est << std::endl;
	  sum += est;
	  n++;
	}
      }
      limit = limit/1.5;
    }
    return (size_t)(sum/n);
  }

  size_t f1() const {
    size_t R = size*32;
    double sum = 0;
    int n = 0;

    double limit = 0.2;
    while (n == 0) {
      std::cout << "R = " << R << std::endl;
      for (size_t i = 0; i < MAX_TABLE; i++) {
	size_t r1 = 0, r0 = 0;
	for (size_t j = 0; j < R; j++) {
	  uint64_t val = getVal(j,i);
	  if (val == 0) {
	    r0++;
	  }
	  if (val == 1) {
	    r1++;
	  }
	}
	
	if ((r0 <= (1-limit)*R) && (r0 >= limit*R)) {
	  sum += (R-1) * (r1/((double) r0)) * pow(2.0,i);
	  n++;
	}
      }
      limit = limit/1.5;
    }
    
    return (size_t) (sum/n);
  }
  
  std::string report() const {
    std::stringstream s;
    /*
    for (size_t i = 0; i < MAX_TABLE; i++) {
      s << i << " " << M[i]<< endl;
      for (size_t j = 0; j < size; j++) {
	s << bitset<64>(table[i*size+j]) << " ";
      }
      s << endl;
    }
    */ 

    s << "F0 = " << F0() << std::endl;
    s << "f1 = " << f1() << std::endl;
    
    return s.str();
  }
private:

  uint64_t getVal(size_t index, size_t w) const {
    size_t wordindex = w*size + (index>>5);
    size_t bitindex = index & 31;
    uint64_t bitmask = maxVal<<(2*bitindex);
    return (table[wordindex] & bitmask) >> (2*bitindex);
  }

  void setVal(size_t index, size_t w, uint64_t val) {
    if (val > maxVal) {
      val = maxVal;
    }
    size_t wordindex = w*size + (index>>5);
    size_t bitindex = index & 31;
    uint64_t bitmask = maxVal<<(2*bitindex);
    table[wordindex] = (((val & maxVal) << (2*bitindex)) & bitmask) | (table[wordindex] & ~bitmask);
  }

  int seed;
  double e;
  uint64_t *table;
  size_t *M;
  size_t size;
  uint64_t mask;
  const size_t MAX_TABLE;
  const uint64_t maxVal; // has to be a power of 2
};


#endif
