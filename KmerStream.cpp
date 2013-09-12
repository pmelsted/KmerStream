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

#include <omp.h>

#include "common.h"
#include "kseq.h"
#include "StreamCounter.hpp"
#include "RepHash.hpp"
#include "lsb.hpp"
#include "StreamCounter.hpp"

//
#include "KmerIterator.hpp"
#include "Kmer.hpp"
#define SEQAN_HAS_ZLIB 1

// seqan libraries
#include <seqan/bam_io.h>


KSEQ_INIT(gzFile, gzread) 

using namespace std;

typedef vector<pair<string,string> > read_t;


struct ProgramOptions {
  size_t k;
  bool verbose;
  bool bam;
  string output;
  vector<string> files;
  double e;
  size_t threads;
  int seed;
  size_t chunksize;
  ProgramOptions() : k(0), verbose(false), bam(false), e(0.01), seed(0), threads(1), chunksize(100000) {}
};

void PrintUsage() {
  cerr << "KmerStream " << VERSION << endl << endl;
  cerr << "Estimates occurrences of k-mers in fastq or fasta files and saves results" << endl << endl;
  cerr << "Usage: KmerStream [options] ... FASTQ files";
  cerr << endl << endl <<
    "-k, --kmer-size=INT      Size of k-mers" << endl <<
    "-q, --quality-cutoff=INT Keep k-mers with bases above quality threshold (default 0)" << endl <<
    "-o, --output=STRING      Filename for output" << endl <<
    "-e, --error-rate=FLOAT   Error rate guaranteed (default value 0.01)" << endl <<
    "-t, --threads=INT        Number of threads to use (default value 1)" << endl <<
    "-s, --seed=INT           Seed value for the randomness (default value 0, use time based randomness)" << endl <<
    "-b, --bam                Input is in BAM format (default false)" << endl <<
    "    --verbose            Print lots of messages during run" << endl << endl
    ;
}

void ParseOptions(int argc, char **argv, ProgramOptions &opt) {
  int verbose_flag = 0;
  int bam_flag = 0;
  const char* opt_string = "k:o:e:s:bt:";
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose_flag, 1},
    {"bam", no_argument, &bam_flag, 'b'},
    {"threads", required_argument, 0, 't'},
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
    case 't':
      opt.threads = atoi(optarg);
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
  if (bam_flag) {
    opt.bam = true;
  }
}


bool CheckOptions(ProgramOptions &opt) {
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

  if (opt.threads <= 0) {
    cerr << "Threads need to be positive" << endl;
    ret = false;
  } else {
    #ifdef _OPENMP
    if (opt.threads > (size_t) omp_get_max_threads()) {
      cerr << "Cant use " << opt.threads << " threads, lowering to " << omp_get_max_threads() << endl;
      opt.threads = omp_get_max_threads();
    }
    #endif
  }

  return ret;
}



template <typename SP>
void RunFastqStream(const ProgramOptions &opt) {
  std::ios_base::sync_with_stdio(false);
  gzFile fp;  
  kseq_t *seq;
  size_t threads = opt.threads;

  vector<SP> sps(threads,SP(opt)); 
  

  // iterate over all reads
  int l = 1;
  size_t readCount;
  size_t chunk = opt.chunksize;
  read_t reads(chunk);
#ifdef _OPENMP
  omp_set_num_threads(threads);  
#endif 
  for (vector<string>::const_iterator it = opt.files.begin(); it != opt.files.end(); ++it) {
    fp = gzopen(it->c_str(), "r");
    seq = kseq_init(fp);
    
    while (l != 0) { // when l is 0 there are no more reads

      readCount = 0;
      reads.clear();
      while ((l = kseq_read(seq) > 0) && readCount < chunk) {  
	readCount++;
	// seq->seq.s is of length seq->seq.l
	// TODO: parallelize this part
	reads.push_back(make_pair(string(seq->seq.s), string(seq->qual.s)));
	//      sp(seq->seq.s, seq->seq.l, seq->qual.s, seq->qual.l);
      }  
      
      // ok, do this in parallel
#pragma omp parallel 
      {
	size_t threadnum = 0;
#ifdef _OPENMP
	threadnum = omp_get_thread_num();
#endif
	
	
#pragma omp for 
	for (size_t j = 0; j < reads.size(); j++) {
	  sps[threadnum](reads[j].first.c_str(), reads[j].first.size(), reads[j].second.c_str(), reads[j].second.size());
	}
	
      }
    }
    kseq_destroy(seq); 
    gzclose(fp); 
  }
  
  for (size_t i = 1; i < sps.size(); i++) {
    sps[0].join(sps[i]);
  }

  SP& sp = sps[0];

  ofstream of;
  of.open(opt.output.c_str(), ios::out);
  of << sp.report();
  of.close();
};

template <typename SP>
void RunBamStream(const ProgramOptions &opt) {
  // TODO: merge this with seqan libs
  std::ios_base::sync_with_stdio(false);
  SP sp(opt);
  // iterate over all reads
  size_t n = 0,t=0;
  for (vector<string>::const_iterator it = opt.files.begin(); it != opt.files.end(); ++it) {
    // open file
    seqan::BamStream bs(it->c_str());
    //    seq = kseq_init(fp);
    seqan::BamAlignmentRecord rec;
    while (!atEnd(bs)) {
      t++;
      if (seqan::readRecord(rec, bs) == 0) {
	sp(seqan::toCString(rec.seq), seqan::length(rec.seq),
	   seqan::toCString(rec.qual), seqan::length(rec.qual));
	n++;
      }
    }  

  }
  cout << "Reads in bam files " << n << " " << t << endl;

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

  
  void operator()(const char* s, size_t l, const char* q, size_t ql) {
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
    sc(val);
  }

  string report() {
    return sc.report();
  }

  void join(const ReadHasher& o) {
    sc.join(o.sc);
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
  if (opt.bam) {
    cout << "running bam file" << endl;
    RunBamStream<ReadHasher>(opt);
  } else {
    RunFastqStream<ReadHasher>(opt);
  }
}


