#include "RepHash.hpp"
#include "mersennetwister.h"

void RepHash::seed(int s) {
    if (s == 0) {
      s = (int) time(NULL);
    }
    // randomize hashvalues
    MTRand mtr(s);
    for (size_t i = 0; i < 32; i++) {
      hvals[i].hi = (((uint64_t)mtr.randInt()) << 32) | mtr.randInt();
      hvals[i].lo = (((uint64_t)mtr.randInt()) << 32) | mtr.randInt();
    }
    // reset state
    h = state_t();
    ht = state_t(); 
  }
