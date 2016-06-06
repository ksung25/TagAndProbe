#include "CEffZFitter.hh"
//#include "CEffPlotter.hh"
#include "KStyle.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings and constants
  //==============================================================================================================

  // handle input arguments
  const std::string conf      = argv[1];         // input configuration file
  const std::string sigpass   = argv[2];	    // signal model for PASS sample
  const std::string bkgpass   = argv[3];	    // background model for PASS sample
  const std::string sigfail   = argv[4];	    // signal model for FAIL sample
  const std::string bkgfail   = argv[5];	    // background model for FAIL sample
  const std::string outdir    = argv[6];        // output directory
  const int         doPU      = atoi(argv[7]);	// PU re-weighting mode
  const int         charge    = atoi(argv[8]);	// probe charge requirement (0, -1, +1)
  const std::string temfname  = argv[9];        // ROOT file for generating MC-based templates
  const int         numtoys   = atoi(argv[10]); // number of toys
  
  // other settings
  const double       massLo    = 60;
  const double       massHi    = 120;
  const double       fitMassLo = massLo;
  const double       fitMassHi = massHi;
  const int          uncMethod = 0;
  const std::string  pufname   = doPU ? "PUWeights_2012.root" : "none";
  const unsigned int runNumLo  = 0;
  const unsigned int runNumHi  = 999999;

  const std::string sigRefDir = "-";
  const std::string bkgRefDir = "-";
  
  std::cout << std::endl;
  std::cout << " <> Using binning configuration: " << conf << std::endl;
  std::cout << std::endl;

  KStyle();  
  CEffZFitter generator;
  generator.initialize(conf, sigpass, bkgpass, sigfail, bkgfail, outdir, temfname, sigRefDir, bkgRefDir, massLo, massHi, fitMassLo, fitMassHi, uncMethod, pufname, charge, runNumLo, runNumHi);
  generator.generateToys(numtoys);

  std::cout << " <> Generated " << numtoys << " toy data files in toys/" << outdir << std::endl;

  return 0;
}
