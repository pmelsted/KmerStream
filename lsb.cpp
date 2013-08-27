#include <stdint.h>
#include "lsb.hpp"

const uint64_t index64[64] = {
   63,  0, 58,  1, 59, 47, 53,  2,
   60, 39, 48, 27, 54, 33, 42,  3,
   61, 51, 37, 40, 49, 18, 28, 20,
   55, 30, 34, 11, 43, 14, 22,  4,
   62, 57, 46, 52, 38, 26, 32, 41,
   50, 36, 17, 19, 29, 10, 13, 21,
   56, 45, 25, 31, 35, 16,  9, 12,
   44, 24, 15,  8, 23,  7,  6,  5
};
 
/**
 * bitScanForward
 * @author Martin Läuter (1997)
 *         Charles E. Leiserson
 *         Harald Prokop
 *         Keith H. Randall
 * "Using de Bruijn Sequences to Index a 1 in a Computer Word"
 * @param bb bitboard to scan
 * @precondition 
 * @return index (0..63) of least significant one bit, returns 63 if b is 0
 */
uint64_t bitScanForward(uint64_t bb) {
   const uint64_t debruijn64 = 0x07EDD5E59A4E28C2ULL;
   return index64[((bb & -bb) * debruijn64) >> 58];
}
