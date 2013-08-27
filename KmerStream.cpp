#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <zlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <fstream>
#include <stdint.h>
#include <sstream>
#include <bitset>

#include "common.h"
#include "kseq.h"
#include "StreamCounter.hpp"
#include "rephash.h"
// #include "cyclichash.h"
#include "lsb.hpp"

KSEQ_INIT(gzFile, gzread) 

using namespace std;


struct ProgramOptions {
  size_t k;
  bool verbose;
  string output;
  vector<string> files;
  double e;
  ProgramOptions() : k(0), verbose(false), e(0.01) {}
};

void PrintUsage() {
  cerr << "KmerStream " << VERSION << endl << endl;
  cerr << "Estimates occurrences of k-mers in fastq or fasta files and saves results" << endl << endl;
  cerr << "Usage: KmerStream [options] ... FASTQ files";
  cerr << endl << endl <<
    "-k, --kmer-size=INT     Size of k-mers" << endl <<
    "-o, --output=STRING     Filename for output" << endl <<
    "-e, --error-rate=FLOAT  Error rate guaranteed (default value 0.01)" << endl <<
    "    --verbose           Print lots of messages during run" << endl << endl
    ;
}

void ParseOptions(int argc, char **argv, ProgramOptions &opt) {
  int verbose_flag = 0;
  const char* opt_string = "k:o:e:";
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose_flag, 1},
    {"kmer-size", required_argument, 0, 'k'},
    {"error-rate", required_argument, 0, 'e'},
    {"output", required_argument, 0, 'o'},
    {0,0,0,0}
  };

  int option_index = 0; 
  int c;

  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 'k':
      opt.k = atoi(optarg);
      break;
    case 'o':
      opt.output = optarg;
      break;
    case 'e':
      opt.e = atof(optarg);
      break;
    default: break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
  
  if (verbose_flag) {
    opt.verbose = true;
  }
}


bool CheckOptions(const ProgramOptions &opt) {
  bool ret = true;


  if (opt.k <= 0) { 
    cerr << "Error, invalid value for kmer-size: " << opt.k << endl;
    cerr << "Value must be at least 1" << endl;
    ret = false;
  }

  if (opt.files.size() == 0) {
    cerr << "Need to specify files for input" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    vector<string>::const_iterator it;
    int intStat;
    for(it = opt.files.begin(); it != opt.files.end(); ++it) {
      intStat = stat(it->c_str(), &stFileInfo);
      if (intStat != 0) {
	cerr << "Error: file not found, " << *it << endl;
	ret = false;
      }
    }
  }
  return ret;
}

template <typename SP>
void RunKmerStream(const ProgramOptions &opt) {
  gzFile fp;  
  kseq_t *seq;  
  SP sp(opt);
  // iterate over all reads
  int l;
  for (vector<string>::const_iterator it = opt.files.begin(); it != opt.files.end(); ++it) {
    fp = gzopen(it->c_str(), "r");
    seq = kseq_init(fp); // STEP 3: initialize seq  
    while ((l = kseq_read(seq)) != -1) {  // don't care about quality strings
      // seq->seq.s is of length seq->seq.l
      // TODO: parallelize this part
      sp(seq->seq.s, seq->seq.l);
    }  

    kseq_destroy(seq); // STEP 5: destroy seq  
    gzclose(fp); // STEP 6: close the file handler  
  }
  

  ofstream of;
  of.open(opt.output.c_str(), ios::out);
  of << sp.report();
  of.close();
}

class StreamCounter{
public:
  StreamCounter(const ProgramOptions &opt) : MAX_TABLE(32), k(opt.k), e(opt.e), hf(opt.k) {
    size = 1<<7; // 1024 bits
    mask = (size<<6) -1; // 2^6=64 bits per block
    M = new size_t[MAX_TABLE];
    memset(M,0,MAX_TABLE);
    table = new uint64_t[size*MAX_TABLE];
    memset(table, 0, size*MAX_TABLE);
    
  }
  ~StreamCounter() { 
    delete[] table;
    delete[] M;
  }
  
  void operator()(const char* s, size_t l) {
    // create hashes for all k-mers
    // operate on hashes
    if (l < k) {
      return;
    } 
    hf.init(s);
    handle(hf.hash());
    for (size_t i = k; i < l; i++) {
      hf.update(s[i-k],s[i]);
      handle(hf.hash());
    }
  }
  
  void handle(uint64_t val) {
    // val is XXX .. XXX1000.. (w-1 times) ..00 0
    size_t w = 1 + bitScanForward(val); // 1-based index of first 1
    if (w >= MAX_TABLE) {
      w = MAX_TABLE-1;
    }
    if (M[w] == size*64) {
      return;
    }
    // val is now XXXX...XX, random
    val = val >> w; // shift away pattern
    
    uint64_t i = val & mask;
    uint64_t index = w*size+(i>>6);
    uint64_t bit = 1ULL<<(i&63);
    if ((table[index] & bit) == 0) {
      table[index] |= bit; // i>>6 is index, 63=2^6-1 is 0..63 shift amount
      M[w]++;
    }
  }

  StreamCounter& join(StreamCounter &o) { return *this;}
  string report() {
    stringstream s;
    for (size_t i = 0; i < MAX_TABLE; i++) {
      s << i << " " << M[i]<< endl;
      for (size_t j = 0; j < size; j++) {
	s << bitset<64>(table[i*size+j]) << " ";
      }
      s << endl;
    }
    return s.str();
  }
private:
  double e;
  size_t k;
  uint64_t *table;
  size_t *M;
  size_t size;
  uint64_t mask;
  RepHash hf;
  const size_t MAX_TABLE;
};





int main(int argc, char** argv) {
  ProgramOptions opt;
  ParseOptions(argc,argv,opt);
  if (!CheckOptions(opt)) {
    PrintUsage();
    exit(1);
  }
  
  RunKmerStream<StreamCounter>(opt);
}


