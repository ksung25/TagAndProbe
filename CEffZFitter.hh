#ifndef CEFFZFITTER_HH
#define CEFFZFITTER_HH

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "ZSignals.hh"
#include "ZBackgrounds.hh"
#include "RooRealVar.h"

class TTree;
class TCanvas;
class TGraphAsymmErrors;
class TH1D;
class TH2D;
class RooFitResult;
class RooAddPdf;
class RooExtendPdf;
class CSignalModel;
class CBackgroundModel;

class CEffZFitter
{
public:
  CEffZFitter();
  ~CEffZFitter();
  
  void initialize(
    const std::string conf,
    const std::string sigpass,
    const std::string bkgpass,
    const std::string sigfail,
    const std::string bkgfail,
    const std::string outdir,
    const std::string temfname,
    const std::string sigRefDir,
    const std::string bkgRefDir,
    const double massLo,
    const double massHi,
    const double fitMassLo,
    const double fitMassHi,
	const int uncMethod,
    const std::string pufname,
    const int charge,
	const unsigned int runNumLo,
    const unsigned int runNumHi
  ); 
  
  void computeEff();

  void generateToys(
    const int numtoys
  );
  void fitToy(
    const std::string methodname,
    const int toynum
  );
  void prepareTrees(const std::string infname);
 
    
protected:

  // parse binning configuration file
  void parseConf(const std::string conf);
  
  // generate MC-based signal templates
  void makeBinnedTemplates(const std::string temfname, const int charge, TH1D *puWeights);
  void makeUnbinnedTemplates(const std::string temfname, const int charge);
  
  // make efficiency graph
  TGraphAsymmErrors* makeEffGraph(const std::vector<double> &edgesv, 
                                  const std::vector<TTree*> &passv,
				  const std::vector<TTree*> &failv,
                                  const std::string name);
  
  // make 2D efficiency map
  void makeEffHist2D(TH2D *hEff, TH2D *hErrl, TH2D *hErrh, 
                     const std::vector<TTree*> &passv,
		     const std::vector<TTree*> &failv,
                     const std::string name);
  
  // efficiency computation
  void performCount(double &resEff, double &resErrl, double &resErrh,
                    const int ibin,
		    const double xbinLo, const double xbinHi,
		    const double ybinLo, const double ybinHi,
                    TTree *passTree, TTree *failTree, 
                    const std::string name, TCanvas *cpass, TCanvas *cfail);
  void performFit(
    double &resEff,
    double &resErrl,
    double &resErrh,
    const int ibin,
	const double xbinLo,
    const double xbinHi,
	const double ybinLo,
    const double ybinHi,
    TTree *passTree,
    TTree *failTree,
    const std::string name,
    TCanvas *cpass,
    TCanvas *cfail,
    unsigned int numThreads=1
  );
  
  // create and destroy PDF models
  void createModels(double ptMin=0, double ptMax=8000, const std::string name="", const int ibin=0);
  void destroyModels();
  // print correlations of fit parameters
  void printCorrelations(std::ostream& os, RooFitResult* res);
  
  // parse fit results file
  void parseFitResults(std::ifstream &ifs, double &eff, double &errl, double &errh);

  // create HTML page
  void makeHTML();
  void makeHTML(const std::string name, const unsigned int nbins);
  
  ///// data members /////
  
  bool fIsInitialized, fPreparedTrees;
  
  // signal and background models
  RooRealVar m = RooRealVar("m", "mass", 60,120);
  CSignalModel::signalType fSigPass = CSignalModel::kNone, fSigFail = CSignalModel::kNone;
  CBackgroundModel::backgroundType fBkgPass = CBackgroundModel::kNone, fBkgFail = CBackgroundModel::kNone;
  CSignalModel     *sigModPass = 0;
  CBackgroundModel *bkgModPass = 0;
  CSignalModel     *sigModFail = 0;
  CBackgroundModel *bkgModFail = 0;
  RooExtendPdf *esignalPass=0, *ebackgroundPass=0, *esignalFail=0, *ebackgroundFail=0;
  RooAddPdf *modelPass=0, *modelFail=0;
  
  // normalization variables for the efficiency extraction
  RooRealVar Nsig     = RooRealVar("Nsig","Signal Yield", 1, 0, 9999999);
  RooRealVar eff      = RooRealVar("eff","Efficiency", 1, 0, 1);
  RooRealVar NbkgPass = RooRealVar("NbkgPass","Background count in PASS sample", 1, 0, 9999999);
  RooRealVar NbkgFail = RooRealVar("NbkgFail","Background count in FAIL sample", 1, 0, 9999999);  
  RooFormulaVar NsigPass = RooFormulaVar("NsigPass","eff*Nsig",RooArgList(eff,Nsig));
  RooFormulaVar NsigFail = RooFormulaVar("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));
  
  double fMassLo, fMassHi;        // signal extraction mass window  
  double fFitMassLo, fFitMassHi;  // fit mass window
  unsigned int fRunNumLo, fRunNumHi;    // run range
  int fCharge;
  
  // efficiency uncertainty calculation method
  // method: 0 -> Clopper-Pearson
  //         1 -> Feldman-Cousins
  int fUncMethod;
  
  // output directory for results
  std::string fOutputDir;
  std::string fFitResDir="";
  std::string fMethodName="";
  std::string fSigRefDir="";
  std::string fBkgRefDir;

  // bin edges for kinematic
  std::vector<double> fPtBinEdgesv;
  std::vector<double> fEtaBinEdgesv;
  std::vector<double> fPhiBinEdgesv;
  std::vector<double> fNPVBinEdgesv;
  std::vector<double> fJetsBinEdgesv;
  std::vector<double> fMETBinEdgesv;
  
  // flags for |eta| and |phi| binning
  bool fDoAbsEta, fDoAbsPhi;
  
  // flags for binnings to compute efficiencies for
  bool fDoPt, fDoEta, fDoPhi, fDoEtaPt, fDoEtaPhi, fDoNPV, fDoJets, fDoMET;
  
  // trees for pass/fail samples
  std::vector<TTree*> fPassTreePtv,     fFailTreePtv;
  std::vector<TTree*> fPassTreeEtav,    fFailTreeEtav;
  std::vector<TTree*> fPassTreePhiv,    fFailTreePhiv;
  std::vector<TTree*> fPassTreeEtaPtv,  fFailTreeEtaPtv;
  std::vector<TTree*> fPassTreeEtaPhiv, fFailTreeEtaPhiv;
  std::vector<TTree*> fPassTreeNPVv,    fFailTreeNPVv;
  std::vector<TTree*> fPassTreeJetsv,   fFailTreeJetsv;
  std::vector<TTree*> fPassTreeMETv,    fFailTreeMETv;

  //weights
  TFile *pufile=0;
  TH1D *puWeights=0;
};

#endif
