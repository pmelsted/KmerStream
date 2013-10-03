#ifndef STREAMCOUNTER_HPP
#define STREAMCOUNTER_HPP

#include <sstream>
#include <stdint.h>
#include <string>
#include <cmath>

#include "lsb.hpp"

size_t roundUpPowerOfTwo(size_t size) {
  size--; 
  size |= size>>1; 
  size |= size>>2; 
  size |= size>>4; 
  size |= size>>8; 
  size |= size>>16; 
  size |= size>>32; 
  size++;
  return size;
}

class StreamCounter{
public:
  StreamCounter(double e_, int seed_) : MAX_TABLE(32), maxVal(15ULL), countWidth(4), countsPerLong(16), e(e_), seed(seed_), sumCount(0) {
    size_t numcounts = (size_t)(48.0/(e*e) + 1); // approx 3 std-dev, true with 0.001 prob.
    F2size  = roundUpPowerOfTwo((size_t)(2.0/(e*e)+1));
    F2table = new uint64_t[F2size];
    memset(F2table,0,F2size * sizeof(uint64_t));

    if (numcounts < 8192) { numcounts = 8192;}
    
    size =(numcounts+countsPerLong-1)/countsPerLong; // size is number of uint64_t's use
    size = roundUpPowerOfTwo(size);
    mask = (size * countsPerLong) -1;
    maxCount = size * countsPerLong * maxVal;

    M = new size_t[MAX_TABLE];
    memset(M,0,MAX_TABLE*sizeof(size_t));

    table = new uint64_t[size*MAX_TABLE];
    memset(table, 0, size*MAX_TABLE*sizeof(uint64_t));
  }

  StreamCounter(const StreamCounter& o) : seed(o.seed), e(o.e), size(o.size), maxCount(o.maxCount), mask(o.mask), MAX_TABLE(o.MAX_TABLE), countWidth(o.countWidth), countsPerLong(o.countsPerLong), maxVal(o.maxVal), sumCount(0), F2size(o.F2size) {
    // copy constructor, creates object of same size with same seed, but empty data
    M = new size_t[MAX_TABLE];
    memset(M,0,MAX_TABLE*sizeof(size_t));

    table = new uint64_t[size*MAX_TABLE];
    memset(table, 0, size*MAX_TABLE*sizeof(uint64_t));

    F2table = new uint64_t[F2size];
    memset(F2table,0,F2size * sizeof(uint64_t));
  }

  ~StreamCounter() { 
    delete[] F2table;
    delete[] table;
    delete[] M;
  }
  
  int getSeed() const {
    return seed;
  }
    
  void operator()(uint64_t hashval) {
    sumCount++;

    // F2 estimation
    ++F2table[(hashval & (F2size-1))];

    // hashval is XXX .. XXX1000.. (w-1 times) ..00 0
    size_t w = bitScanForward(hashval); // 1-based index of first 1

    if (w >= MAX_TABLE) {
      w = MAX_TABLE-1;
    }

    if (M[w] == size*countsPerLong*maxVal) {
      return;
    }

    // hashval is now XXXX...XX, random
    uint64_t hval = hashval >> (w+1); // shift away pattern of XXX10000
    
    uint64_t index = hval & mask;
    uint64_t val = getVal(index,w);
    if (val != maxVal) { // max count
      setVal(index,w,val+1);
      M[w]++;
    }
  }

  bool join(const StreamCounter &o) {
    if (size != o.size || F2size != o.F2size || seed != o.seed) {
      return false;
    }

    for (size_t i = 0; i < MAX_TABLE; i++) {
      for (size_t j = 0; j < size*countsPerLong; j++) {
	setVal(j,i,getVal(j,i)+o.getVal(j,i));
      }
    }

    sumCount += o.sumCount;
    for (size_t i = 0; i < F2size; i++) {
      F2table[i] += o.F2table[i];
    }

    return true;
  }
  
  size_t F0() const {
    size_t R = size*countsPerLong;
    double sum = 0;

    int n = 0;
    double limit = 0.2;
    while (n == 0 && limit > 1e-8) {
      for (size_t i = 0; i < MAX_TABLE; i++) {
	size_t ts = 0;
	for (size_t j = 0; j < R; j++) {
	  if (getVal(j,i) > 0) {
	    ts++;
	  }
	}
	
	if (ts <= (1-limit)*R && ts >= limit*R) {
	  double est = (log(1.0-ts/((double) R))/log(1.0-1.0/R)) * pow(2.0,i+1);
	  //std::cout << i << " " << ts << " " << limit<<  " " << est << std::endl;
	  sum += est;
	  n++;
	  break;
	}
      }
      limit = limit/1.5;
    }
    return (size_t)(sum/n);
  }

  size_t f1() const {
    size_t R = size*countsPerLong;
    double sum = 0;
    int n = 0;

    double limit = 0.2;
    while (n == 0 && limit > 1e-8) {
      //std::cout << "R = " << R << std::endl;
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
	  sum += (R-1) * (r1/((double) r0)) * pow(2.0,i+1);
	  n++;
	  break;
	}
      }
      limit = limit/1.5;
    }
    
    return (size_t) (sum/n);
  }

  /*
 size_t f2() const {
    size_t R = size*countsPerLong;
    double sum = 0;
    int n = 0;

    double limit = 0.2;
    while (n == 0 && limit > 1e-8) {
      //std::cout << "R = " << R << std::endl;
      for (size_t i = 0; i < MAX_TABLE; i++) {
	size_t r2 = 0, r1 = 0, r0 = 0;
	for (size_t j = 0; j < R; j++) {
	  uint64_t val = getVal(j,i);
	  if (val == 0) {
	    r0++;
	  }
	  if (val == 1) {
	    r1++;
	  }
	  if (val == 2) {
	    r2++;
	  }
	}
	
	if ((r0 <= (1-limit)*R) && (r0 >= limit*R)) {
	  double x1 = (r1/((double) r0));
	  double x2 = (r2/((double) r0) - x1*x1/2.0 );
	  std::cout << "r0,r1,r2 = " << r0 << ", " << r1 << ", " << r2 << std::endl;
	  sum += x2*(R-1)*pow(2.0,i+1);
	  n++;
	  break;
	}
      }
      limit = limit/1.5;
    }
    
    return (size_t) (sum/n);
  }
 size_t f3() const {
    size_t R = size*countsPerLong;
    double sum = 0;
    int n = 0;

    double limit = 0.2;
    while (n == 0 && limit > 1e-8) {
      //std::cout << "R = " << R << std::endl;
      for (size_t i = 0; i < MAX_TABLE; i++) {
	size_t r3 = 0, r2 = 0, r1 = 0, r0 = 0;
	for (size_t j = 0; j < R; j++) {
	  uint64_t val = getVal(j,i);
	  if (val == 0) {
	    r0++;
	  }
	  if (val == 1) {
	    r1++;
	  }
	  if (val == 2) {
	    r2++;
	  }
	  if (val == 3) {
	    r3++;
	  }
	}
	
	if ((r0 <= (1-limit)*R) && (r0 >= limit*R)) {
	  double x1 = (r1/((double) r0));
	  double x2 = (r2/((double) r0) - x1*x1/2.0 );
	  double x3 = (r3/((double) r0) - x1*x1*x1/6.0 - x1*x1*x2/2.0);
	  std::cout << "r0,r1,r2 = " << r0 << ", " << r1 << ", " << r2 << ", " << r3 << std::endl;
	  sum += x3*(R-1)*pow(2.0,i+1);
	  n++;
	  break;
	}
      }
      limit = limit/1.5;
    }
    
    return (size_t) (sum/n);
    }*/

  
  std::string report() const {
    std::stringstream s;
    /*
    size_t R = size*countsPerLong;
    s << "R = " << R << std::endl;
    for (size_t i = 0; i < MAX_TABLE; i++) {
      int count[4] = {0,0,0,0};
      s << "M[" << i << "] =  " << M[i]<< std::endl;
      for (size_t j = 0; j < R; j++) {
	uint64_t val = getVal(j,i);
	count[val]++;
      }
      s << count[0] << ", " << count[1] << ", " << count[2] << ", " << count[3] << std::endl;
    }
    s << std::endl << F0() << std::endl;
    */
    
    s << "F0 = " << F0() << std::endl;
    s << "f1 = " << f1() << std::endl;
  //s << "f2 = " << f2() << std::endl;
  //s << "f3 = " << f3() << std::endl;
    s << "F1 = " << sumCount << std::endl;
    s << "F2 = " << F2() << std::endl;
    
    
    return s.str();
  }

  size_t F2() const {
    // based on Thorup and Zhangs, F2 estimator
    double sum=0,sqsum=0;
    for (size_t i = 0; i < F2size; i++) {
      double c = (double) F2table[i];
      sum += c;
      sqsum += c*c;
    }
    return (size_t)(sqsum + (sqsum - sum*sum)/F2size);
  }

private:

  uint64_t getVal(size_t index, size_t w) const {
    size_t wordindex = w*size + (index/countsPerLong);
    size_t bitindex = index & (countsPerLong-1) ;
    uint64_t bitmask = maxVal<<(countWidth*bitindex);
    return (table[wordindex] & bitmask) >> (countWidth*bitindex);
  }

  void setVal(size_t index, size_t w, uint64_t val) {
    if (val > maxVal) {
      val = maxVal;
    }
    size_t wordindex = w*size + (index/countsPerLong);
    size_t bitindex = index & (countsPerLong-1);
    uint64_t bitmask = maxVal<<(countWidth*bitindex);
    table[wordindex] = (((val & maxVal) << (countWidth*bitindex)) & bitmask) | (table[wordindex] & ~bitmask);
  }

  int seed;
  double e;
  uint64_t *table;
  uint64_t *F2table;
  size_t F2size;
  size_t sumCount;
  size_t *M;
  size_t size;
  size_t maxCount;
  uint64_t mask;
  const size_t MAX_TABLE;
  const size_t countWidth; // number of bits per count, even number
  const size_t countsPerLong;  // fix this
  const uint64_t maxVal; // has to be a power of 2-1
};


#endif
