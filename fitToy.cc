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
  const std::string methodname= argv[1];
  const std::string conf      = argv[2];         // input configuration file
  const std::string sigpass   = argv[3];	     // signal model for PASS sample
  const std::string bkgpass   = argv[4];	     // background model for PASS sample
  const std::string sigfail   = argv[5];	     // signal model for FAIL sample
  const std::string bkgfail   = argv[6];	     // background model for FAIL sample
  const std::string outdir    = argv[7];         // output directory
  const int         doPU      = atoi(argv[8]);	// PU re-weighting mode
  const int         charge    = atoi(argv[9]);	// probe charge requirement (0, -1, +1)
  const std::string temfname  = argv[10];        // ROOT file for generating MC-based templates
  const std::string sigRefDir = argv[11];        // reference dir for fixed signal parameters
  const std::string bkgRefDir = argv[12];        // reference dir for fixed background parameters
  const int         toynum    = atoi(argv[13]); // toy number to fit
  
  // other settings
  const double       massLo    = 60;
  const double       massHi    = 120;
  const double       fitMassLo = massLo;
  const double       fitMassHi = massHi;
  const int          uncMethod = 0;
  const std::string  pufname   = doPU ? "PUWeights_2012.root" : "none";
  const unsigned int runNumLo  = 0;
  const unsigned int runNumHi  = 999999;
  
  std::cout << std::endl;
  std::cout << " <> Using binning configuration: " << conf << std::endl;
  std::cout << std::endl;

  KStyle();  
  CEffZFitter fitter;
  fitter.initialize(conf, sigpass, bkgpass, sigfail, bkgfail, outdir, temfname, sigRefDir, bkgRefDir, massLo, massHi, fitMassLo, fitMassHi, uncMethod, pufname, charge, runNumLo, runNumHi);
  if(toynum>=0) {
    std::cout << " <> Fitting toys in toys/" << outdir << std::endl;
    fitter.fitToy(methodname, toynum);
  } else {
    std::cout << " <> Initialized without fitting any toys " << outdir << std::endl;
  }

  return 0;
}
