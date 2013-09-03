#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>

#include <sstream>
#include <vector>
#include <string>

#include <stdint.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include "Common.hpp"


#include "fastq.hpp"
#include "Kmer.hpp"
#include "lsb.hpp"


// structs for getopt

struct CountStream_ProgramOptions {
  size_t k;
  size_t m;
  string output;
  bool verbose;
  vector<string> files;
  CountStream_ProgramOptions() : k(0), m(1024), verbose(false) {}
};

void CountStream_PrintUsage() {
  cerr << "StreamCounter " << SC_VERSION << endl << endl;
  cerr << "Estimates number of k-mers in sequence files" << endl << endl;
  cerr << "Usage: StreamCounter stream [options] ... FASTQ files";
  cerr << endl << endl <<
    "-k, --kmer-size=INT             Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl << 
    "-m, --bucket-size=INT           Number of buckets, should be a power of 2, recommended values" << endl <<
    "                                are 256 (7% error), 1024 (2% error, default) and 65536 (.5% error)" << endl <<
    "-o, --output=STRING             Filename for output" << endl <<
    "    --verbose                   Print lots of messages during run" << endl << endl
    ;
}




void CountStream_ParseOptions(int argc, char **argv, CountStream_ProgramOptions &opt) {
  int verbose_flag = 0;
  const char* opt_string = "k:m:o:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
    {"kmer-size", required_argument, 0, 'k'},
    {"bucket-size", required_argument, 0, 'm'},
    {"output", required_argument, 0, 'o'},
    {0,0,0,0}
  };

  int option_index = 0;
  int c;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);
    //cout << "debug: c="<<c<<", optarg="<<optarg << endl;
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
    case 'm':
      opt.m = atoi(optarg);
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


bool CountStream_CheckOptions(CountStream_ProgramOptions &opt) {
  bool ret = true;

  if (opt.k <= 0 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error, invalid value for kmer-size: " << opt.k << endl;
    cerr << "Values must be between 1 and " << (MAX_KMER_SIZE-1) << endl;
    ret = false;
  }

  if (opt.m <= 0) {
    cerr << "Error, invalid value for bucket-size: " << opt.m << endl;
    cerr << "Values must be positive integers and powers of two" << endl;
    ret = false;
  }
  if (opt.m < 16) {
    opt.m = 16;
  }

  if (opt.m > (1<<16)) {
    opt.m = (1<<16);
  }
  
  if ((opt.m & (opt.m-1)) != 0) {
    cerr << "Error, m should be a power of 2: " << opt.m << endl;
    size_t x = opt.m-1;
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x | (x >> 16);
    opt.m = x+1;
    cerr << "value modified to m = " << opt.m << endl;
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
  
 
  //TODO: check if we have permission to write to outputfile
  
  return ret;

}

void CountStream_PrintSummary(const CountStream_ProgramOptions &opt) {
 
}


void CountStream_Normal(const CountStream_ProgramOptions &opt) {
  // create hash table and bloom filter
 
  char name[8196],s[8196];
  size_t name_len,len;

  size_t k = opt.k;
  size_t m = opt.m;
  size_t phi_inv = 20; // TODO make this optional
  size_t b = 0;
  while ((m & (1<<b)) == 0) {
    b++;
  }
  cerr << "m,b = " << m << "," << b << endl;

  uint32_t seed = (uint32_t) time(NULL);


  uint64_t n_read = 0;
  uint64_t num_kmers = 0;  

  size_t table_size = phi_inv * m;

  uint8_t *V = new uint8_t[table_size*64];
  memset(V,0,table_size*64);
  int *M = new int[m];
  for (size_t i = 0; i < m; i++) {
    M[i] = 0;
  }

  // loops over all files
  FastqFile FQ(opt.files);

  // for each read
  while (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0) {
    // TODO: add code to handle N's, currently all N's are mapped to A
    Kmer km(s);
    Kmer tw = km.twin();
    for (size_t i = 0; i <= len-k; ++i) {
      ++num_kmers;
      if (i > 0) {
	km = km.forwardBase(s[i+k-1]);
	tw = tw.backwardBase(s[i+k-1]);
      }
      Kmer rep = (km < tw) ? km : tw;
      // process rep
      size_t h = rep.hash(seed);
      size_t j = h >> (64-b);
      size_t w = 1+bitScanForward(h); // Hyperloglog needs 1-based index, TODO shift h prior

      if (w > 64-b) {
	w = 64-b;
      }
      if (M[j] < w) {
	M[j] = w;
      }

      // TODO: figure out how to put in correct bin
      /*      if (V[(64*j+w)] < 2) {
	V[(64*j+w)]++;
	}*/

    }
    ++n_read;

    if (opt.verbose && n_read % 1000000 == 0) {
      cerr << "processed " << n_read << " reads" << endl;
    }
  }

  if (opt.verbose) {
    cerr << "processed " << num_kmers << " kmers in " << n_read  << " reads"<< endl;
  }

  
  /*int *M = new int[m];
  for (size_t i = 0; i < m; i++) {
    M[i] = 0;
    for (size_t j = 0; j < 64; j++) {
      if (V[(64*i+j)] > 0) {
	M[i] = (int)j;
      }
    }
  }
  */

  double Zinv = 0.0;
  for (size_t i = 0; i < m; i++) {
    Zinv += pow(2.0,-M[i]);
  }
  
  // compute a_m
  double a_m = 0.0;
  if (m == 16) {
    a_m = 0.673;
  } else if (m == 32) {
    a_m = 0.697;
  } else if (m == 64) {
    a_m = 0.709;
  } else {
    a_m = 0.7213/(1+1.079/m);
  }

  double E = (a_m*m)*m/Zinv;
  cout << "a_m      " << a_m << endl;
  cout << "Estimate " << E << endl;
  cout << "(int)    " << ((size_t) E) << endl;

  // cleanup
  //delete[] V;
  delete[] M;
  
}

void CountStream(int argc, char **argv) {
  
  CountStream_ProgramOptions opt;
  CountStream_ParseOptions(argc,argv,opt);

  if (argc < 2) {
    CountStream_PrintUsage();
    exit(1);
  }
  
  if (!CountStream_CheckOptions(opt)) {
    CountStream_PrintUsage();
    exit(1);
  }
  
  // set static global k-value
  Kmer::set_k(opt.k);

  if (opt.verbose) {
    CountStream_PrintSummary(opt);
  }

  CountStream_Normal(opt);

}
