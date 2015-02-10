#include <iostream>
#include <string>
#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <algorithm>

#include "rephash.h"

int main(int argc, char **argv) {
  size_t k = atoi(argv[1]);

  const char *fw = "TGATCTTTCAGATTGTAGAGTTTCATTTAGTTTACCAGTACTCGTGCGCCCGCCGAATCCAGGCGTCAAA";
  const char *bw = "TTTGACGCCTGGATTCGGCGGGCGCACGAGTACTGGTAAACTAAATGAAACTCTACAATCTGAAAGATCA";
  size_t l = std::string(fw).size();
  std::vector<uint64_t> fh,bh;

  RepHash hf(k);
  hf.seed(42);
  hf.init(fw);
  fh.push_back(hf.hash());

  for (size_t i = k; i < l; i++) {
    hf.update(fw[i-k],fw[i]);
    fh.push_back(hf.hash());
  }

  hf.init(bw);
  bh.push_back(hf.hash());

  for (size_t i = k; i < l; i++) {
    hf.update(bw[i-k],bw[i]);
    bh.push_back(hf.hash());
  }
  reverse(bh.begin(),bh.end());

  for (size_t i = 0; i < fh.size(); i++) {
    std::cout << fh[i] << std::endl << bh[i] << std::endl;
  }


  std::cout << (fh == bh) << std::endl;


}
