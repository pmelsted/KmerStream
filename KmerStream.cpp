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
#include "lsb.hpp"
#include "StreamCounter.hpp"

//
#include "KmerIterator.hpp"
#include "Kmer.hpp"


KSEQ_INIT(gzFile, gzread) 

using namespace std;


struct ProgramOptions {
  size_t k;
  bool verbose;
  string output;
  vector<string> files;
  double e;
  int seed;
  ProgramOptions() : k(0), verbose(false), e(0.01), seed(0) {}
};

void PrintUsage() {
  cerr << "KmerStream " << VERSION << endl << endl;
  cerr << "Estimates occurrences of k-mers in fastq or fasta files and saves results" << endl << endl;
  cerr << "Usage: KmerStream [options] ... FASTQ files";
  cerr << endl << endl <<
    "-k, --kmer-size=INT     Size of k-mers" << endl <<
    "-o, --output=STRING     Filename for output" << endl <<
    "-e, --error-rate=FLOAT  Error rate guaranteed (default value 0.01)" << endl <<
    "-s, --seed=INT          Seed value for the randomness (default value 0, use time based randomness)" << endl <<
    "    --verbose           Print lots of messages during run" << endl << endl
    ;
}

void ParseOptions(int argc, char **argv, ProgramOptions &opt) {
  int verbose_flag = 0;
  const char* opt_string = "k:o:e:s:";
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose_flag, 1},
    {"kmer-size", required_argument, 0, 'k'},
    {"error-rate", required_argument, 0, 'e'},
    {"seed", required_argument, 0, 's'},
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
    case 's':
      opt.seed = atoi(optarg);
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
void RunStream(const ProgramOptions &opt) {
  std::ios_base::sync_with_stdio(false);
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
};


class ReadHasher {
public:
  ReadHasher(const ProgramOptions &opt) : k(opt.k), hf(opt.k), sc(opt.e, opt.seed) {
    if (opt.seed != 0) {
      hf.seed(opt.seed);
    }
  }

  
  void operator()(const char* s, size_t l) {
    // create hashes for all k-mers
    // operate on hashes
    /*
    if (l < k) {
      return;
    } 
    hf.init(s);
    handle(hf.hash());
    for (size_t i = k; i < l; i++) {
      hf.update(s[i-k],s[i]);
      handle(hf.hash());
      }*/
    KmerIterator it(s),it_end;
    for (; it != it_end; ++it) {
      handle(it->first.rep().hash());
    }
    
  }
  
  void handle(uint64_t val) {
    //sc(val);
    cout << val << "\n";
  }

  string report() {
    return sc.report();
  }

private:
  RepHash hf;
  size_t k;
  StreamCounter sc;
  //F2Counter f2;
  //SumCounter sum;
};





int main(int argc, char** argv) {
  ProgramOptions opt;
  ParseOptions(argc,argv,opt);
  if (!CheckOptions(opt)) {
    PrintUsage();
    exit(1);
  }
  
  //
  Kmer::set_k(opt.k);
  RunStream<ReadHasher>(opt);
}


