#ifndef REPHASH_H
#define REPHASH_H

#include <stdint.h>
#include <cassert>
#include <time.h>

static const unsigned char twin[32] = {
  0, 20, 2, 7, 4, 5, 6, 3,  
  8,  9, 10, 11, 12, 13, 14, 15, 
  16, 17, 18, 19, 1, 21, 22, 23, 
  24, 25, 26, 27, 28, 29, 30, 31
};

struct state_t {
  uint64_t hi;
  uint64_t lo;
  state_t& operator^=(state_t& o) {
    hi ^= o.hi;
    lo ^= o.lo;
    return *this;
  }
  
  state_t() : hi(0), lo(0) {}
  
};

class RepHash {
 public:

 

	RepHash(int _k) {
		assert(_k>0);
		assert(_k<=63);
		init(_k);
		seed(0);
	}

	RepHash() {
		init(0);
		seed(0);
	}

	void init(int _k) {
		k = _k;
		last1mask =  (1ULL << 63);
		lastkmask =  ((1ULL<<k)-1) << (64-k);
		firstkmask = (1ULL<<k)-1;
		charmask = 31;
	}
  
  void seed(int s);
  
  

  inline void fastleftshiftk(state_t& x) {
    uint64_t upper =    (x.hi & (lastkmask)); // upper k of 64 bits of hi
    x.hi = (x.hi<<k) | ((x.lo & (lastkmask)) >> (64-k));// 
    x.lo = (x.lo<<k) | ((upper)              >> (64-k));
  }

  inline void fastrightshiftk(state_t& x) {
    uint64_t lower =    (x.hi & (firstkmask)); // lower k bits
    x.hi = (x.hi>>k) | ((x.lo & (firstkmask)) << (64-k));
    x.lo = (x.lo>>k) | ((lower)               << (64-k));
  }

  inline void fastleftshift1(state_t& x) {
    uint64_t last1 = (x.hi & last1mask); // last bit of hi
    x.hi = (x.hi<<1) | ((x.lo & last1mask)>>63);
    x.lo = (x.lo<<1) | (last1 >> 63);
  }
  
  inline void fastrightshift1(state_t& x) {
    uint64_t first1 = (x.hi & 1ULL);  // first bit of hi
    x.hi = (x.hi>>1) | ((x.lo & 1ULL) <<63);
    x.lo = (x.lo>>1) | (first1 << 63);
  }

  uint64_t hash() const {
    return h.lo ^ ht.lo;
  }

  void init(const char* _s) {
    h = state_t();
    ht = state_t();
    const unsigned char* s = (const unsigned char*) _s;
    for (size_t i = 0; i < k; i++) {
      fastleftshift1(h);
      state_t hval = hvals[s[i] & charmask];
      h ^= hval;

      fastleftshift1(ht);
      state_t hvalt = hvals[twin[s[k-1-i] & charmask]];
      ht ^= hvalt;
    }
  }

  inline void update(const unsigned char out, const unsigned char in) {
    state_t z(hvals[out & charmask]);
    fastleftshiftk(z);
    fastleftshift1(h);
    state_t hval(hvals[in & charmask]);
    h ^= z;
    h ^= hval;
    
    state_t zt(hvals[twin[in & charmask]]);
    fastleftshiftk(zt);
    state_t hvalt(hvals[twin[out & charmask]]);
    ht ^= hvalt;
    ht ^= zt;
    fastrightshift1(ht);
  }

 private:
  size_t k;
  uint64_t last1mask, lastkmask, firstkmask;
  unsigned char charmask;
  state_t h,ht;
  state_t hvals[32];
};




#endif // REPHASH_H
