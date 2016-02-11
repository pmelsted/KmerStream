#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <sys/stat.h>
#include <fstream>
#include <stdint.h>
#include <sstream>
#include <bitset>
#include <cstdlib>


#include "common.h"
#include "StreamCounter.hpp"


using namespace std;

struct ProgramOptions {
  //size_t k;
  std::string output;
  std::vector<string> files;
  bool verbose;
  ProgramOptions() : verbose(false) {}
};

void PrintUsage() {
  cerr << "StreamJoin " << VERSION << endl << endl;
  cerr << "Creates union of many stream estimates" << endl << endl;
  cerr << "Usage: StreamJoin [options] ... files";
  cerr << endl << endl <<
       "-o, --output=STRING      Filename for output" << endl <<
       "    --verbose            Print output at the end" << endl <<
       endl;

}

void ParseOptions(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;

  const char *opt_string = "o:";
  static struct option long_options[] = {
    {"verbose", no_argument, &verbose_flag, 1},
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
    case 'o':
      opt.output = optarg;
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


bool CheckOptions(ProgramOptions& opt) {
  bool ret = true;

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

  if (opt.output.empty()) {
    if (opt.files.size() != 1) {
      cerr << "Missing output file" << endl;
      ret = false;
    }
  }

  return ret;
}





int main(int argc, char **argv) {
  ProgramOptions opt;
  ParseOptions(argc,argv,opt);
  if (!CheckOptions(opt)) {
    PrintUsage();
    exit(1);
  }

  if (opt.output.empty()) {
    std::string fn = opt.files[0];
    StreamCounter sc(0.01,0);
    sc.loadBinary(fn);
    std::cout << sc.report();
  } else {
    StreamCounter sc(0.01,0);
    StreamCounter join(0.01,0);
    join.loadBinary(opt.files[0]);
    for(size_t i = 1; i < opt.files.size(); i++) {
      sc.loadBinary(opt.files[i]);
      join.join(sc);
    }
    join.writeBinary(opt.output);
    if (opt.verbose) {
      std::cout << join.report();
    }

  }
}
