#include "CEffZFitter.hh"
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TColor.h>
#include <TROOT.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>

#include <cassert>
#include <sstream>
#include <iomanip>
#include <unistd.h>

#include "CPlot.hh"
#include "KStyle.hh"
//#include "ZSignals.hh"
//#include "ZBackgrounds.hh"
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"

// RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooMCStudy.h"

// bin size constants
#define BIN_SIZE_PASS 2
#define BIN_SIZE_FAIL 2

//--------------------------------------------------------------------------------------------------
CEffZFitter::CEffZFitter():
fIsInitialized(false),
fPreparedTrees(false),
//fSigPass      (CSignalModel::kNone),
//fBkgPass      (CBackgroundModel::kNone),
//fSigFail      (CSignalModel::kNone),
//fBkgFail      (CBackgroundModel::kNone),
fMassLo       (60),
fMassHi       (120),
fFitMassLo    (60),
fFitMassHi    (120),
fRunNumLo     (0),
fRunNumHi     (99999999),
fCharge       (1),
fUncMethod    (0),
fOutputDir    ("."),
fDoAbsEta     (false),
fDoAbsPhi     (false),
fDoPt         (false),
fDoEta        (false),
fDoPhi        (false),
fDoEtaPt      (false),
fDoEtaPhi     (false),
fDoNPV        (false),
fDoJets       (false),
fDoMET        (false)
{}

//--------------------------------------------------------------------------------------------------
CEffZFitter::~CEffZFitter()
{
  for(unsigned int i=0; i<fPassTreePtv.size(); i++)     { delete fPassTreePtv[i];     fPassTreePtv[i]=0;  }
  for(unsigned int i=0; i<fPassTreeEtav.size(); i++)    { delete fPassTreeEtav[i];    fPassTreeEtav[i]=0; }
  for(unsigned int i=0; i<fPassTreePhiv.size(); i++)    { delete fPassTreePhiv[i];    fPassTreePhiv[i]=0; }
  for(unsigned int i=0; i<fPassTreeEtaPtv.size(); i++)  { delete fPassTreeEtaPtv[i];  fPassTreeEtaPtv[i]=0; }
  for(unsigned int i=0; i<fPassTreeEtaPhiv.size(); i++) { delete fPassTreeEtaPhiv[i]; fPassTreeEtaPhiv[i]=0; }
  for(unsigned int i=0; i<fPassTreeNPVv.size(); i++)    { delete fPassTreeNPVv[i];    fPassTreeNPVv[i]=0; }
  for(unsigned int i=0; i<fPassTreeJetsv.size(); i++)   { delete fPassTreeJetsv[i];   fPassTreeJetsv[i]=0; }
  for(unsigned int i=0; i<fPassTreeMETv.size(); i++)    { delete fPassTreeMETv[i];    fPassTreeMETv[i]=0; }

  for(unsigned int i=0; i<fFailTreePtv.size(); i++)     { delete fFailTreePtv[i];     fFailTreePtv[i]=0; }
  for(unsigned int i=0; i<fFailTreeEtav.size(); i++)    { delete fFailTreeEtav[i];    fFailTreeEtav[i]=0; }
  for(unsigned int i=0; i<fFailTreePhiv.size(); i++)    { delete fFailTreePhiv[i];    fFailTreePhiv[i]=0; }
  for(unsigned int i=0; i<fFailTreeEtaPtv.size(); i++)  { delete fFailTreeEtaPtv[i];  fFailTreeEtaPtv[i]=0; }
  for(unsigned int i=0; i<fFailTreeEtaPhiv.size(); i++) { delete fFailTreeEtaPhiv[i]; fFailTreeEtaPhiv[i]=0; }
  for(unsigned int i=0; i<fFailTreeNPVv.size(); i++)    { delete fFailTreeNPVv[i];    fFailTreeNPVv[i]=0; }
  for(unsigned int i=0; i<fFailTreeJetsv.size(); i++)   { delete fFailTreeJetsv[i];   fFailTreeJetsv[i]=0; }
  for(unsigned int i=0; i<fFailTreeMETv.size(); i++)    { delete fFailTreeMETv[i];    fFailTreeMETv[i]=0; }
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::initialize(
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
) {
  std::cout << "   [CEffZFitter] Initializing... " << std::endl;

  //parse binning configuration file
  parseConf(conf);
  
  if(sigpass=="BreitWignerConvCrystalBall" ) fSigPass=CSignalModel::kBreitWignerConvCrystalBall ;
  else if(sigpass=="MCTemplateConvGaussian"     ) fSigPass=CSignalModel::kMCTemplateConvGaussian     ;
  else if(sigpass=="MCTemplateConvGaussianInit" ) fSigPass=CSignalModel::kMCTemplateConvGaussianInit ;
  else if(sigpass=="VoigtianCBShape"            ) fSigPass=CSignalModel::kVoigtianCBShape            ;
  else if(sigpass=="MCDatasetConvGaussian"      ) fSigPass=CSignalModel::kMCDatasetConvGaussian      ;
  else if(sigpass=="BWCBPlusVoigt"              ) fSigPass=CSignalModel::kBWCBPlusVoigt              ;
  else if(sigpass=="BWCBPlusVoigtBounded"       ) fSigPass=CSignalModel::kBWCBPlusVoigtBounded       ;
  else fSigPass=CSignalModel::kNone;
  if(sigfail=="BreitWignerConvCrystalBall" ) fSigFail=CSignalModel::kBreitWignerConvCrystalBall ;
  else if(sigfail=="MCTemplateConvGaussian"     ) fSigFail=CSignalModel::kMCTemplateConvGaussian     ;
  else if(sigfail=="MCTemplateConvGaussianInit" ) fSigFail=CSignalModel::kMCTemplateConvGaussianInit ;
  else if(sigfail=="VoigtianCBShape"            ) fSigFail=CSignalModel::kVoigtianCBShape            ;
  else if(sigfail=="MCDatasetConvGaussian"      ) fSigFail=CSignalModel::kMCDatasetConvGaussian      ;
  else if(sigfail=="BWCBPlusVoigt"              ) fSigFail=CSignalModel::kBWCBPlusVoigt              ;
  else if(sigfail=="BWCBPlusVoigtBounded"       ) fSigFail=CSignalModel::kBWCBPlusVoigtBounded       ;
  else fSigFail=CSignalModel::kNone;
  if(bkgpass=="Exponential"   ) fBkgPass=CBackgroundModel::kExponential  ;
  else if(bkgpass=="ErfcExpo"      ) fBkgPass=CBackgroundModel::kErfcExpo     ;
  else if(bkgpass=="ErfcExpoFixed" ) fBkgPass=CBackgroundModel::kErfcExpoFixed;
  else if(bkgpass=="DoubleExp"     ) fBkgPass=CBackgroundModel::kDoubleExp    ;
  else if(bkgpass=="LinearExp"     ) fBkgPass=CBackgroundModel::kLinearExp    ;
  else if(bkgpass=="QuadraticExp"  ) fBkgPass=CBackgroundModel::kQuadraticExp ;
  else fBkgPass=CBackgroundModel::kNone;
  if(bkgfail=="Exponential"   ) fBkgFail=CBackgroundModel::kExponential  ;
  else if(bkgfail=="ErfcExpo"      ) fBkgFail=CBackgroundModel::kErfcExpo     ;
  else if(bkgfail=="ErfcExpoFixed" ) fBkgFail=CBackgroundModel::kErfcExpoFixed;
  else if(bkgfail=="DoubleExp"     ) fBkgFail=CBackgroundModel::kDoubleExp    ;
  else if(bkgfail=="LinearExp"     ) fBkgFail=CBackgroundModel::kLinearExp    ;
  else if(bkgfail=="QuadraticExp"  ) fBkgFail=CBackgroundModel::kQuadraticExp ;
  else fBkgFail=CBackgroundModel::kNone;
  fMassLo    = massLo;
  fMassHi    = massHi;
  fFitMassLo = fitMassLo;
  fFitMassHi = fitMassHi;
  fRunNumLo  = runNumLo;
  fRunNumHi  = runNumHi;
  fCharge    = charge;
  m.setRange(fFitMassLo, fFitMassHi); 
  m.setBins(10000);

  // set up output directory
  fOutputDir = outdir;
  fSigRefDir = sigRefDir;
  fBkgRefDir = bkgRefDir;
  gSystem->mkdir(fOutputDir.c_str(),true);
  CPlot::sOutDir = TString(outdir.c_str()) + TString("/plots");
      
  if(pufname.compare("none")!=0) {
    pufile = new TFile(pufname.c_str());      assert(pufile);
    puWeights = (TH1D*)pufile->Get("pileup"); assert(puWeights); 
  }
  
  // Generate templates from MC if necessary
  // Check if templates file exists and try to make them if they do not exist
  // This will fail if a TNP job is submitted to grid and the template files were not already made
  // interactively and sent with the job.
  // This is necessary because reading the skim file from which the templates are made is I/O
  // intensive across mounted file systems.
  if(
    fSigPass==CSignalModel::kMCTemplateConvGaussian || 
    fSigPass==CSignalModel::kMCTemplateConvGaussianInit ||
    fSigFail==CSignalModel::kMCTemplateConvGaussian || 
    fSigFail==CSignalModel::kMCTemplateConvGaussianInit
  ) {
    string outfile_name = "templates/" + fOutputDir + "/binnedTemplates.root";
    gSystem->mkdir(("templates/"+fOutputDir).c_str(),true);
    int res = access(outfile_name.c_str(), R_OK);
    if(res<0) makeBinnedTemplates(temfname, fCharge, puWeights);
  } else if(
    fSigPass==CSignalModel::kMCDatasetConvGaussian ||
    fSigFail==CSignalModel::kMCDatasetConvGaussian
  ) {
    string outfile_name = "templates/" + fOutputDir + "/unbinnedTemplates.root";
    gSystem->mkdir(("templates/"+fOutputDir).c_str(),true);
    int res = access(outfile_name.c_str(), R_OK);
    if(res<0) makeUnbinnedTemplates(temfname, fCharge);
  }
  fIsInitialized = true;
}
//--------------------------------------------------------------------------------------------------
void CEffZFitter::prepareTrees(
  const std::string infname
)
{
  assert(fIsInitialized);
  std::cout << "   [CEffZFitter] Preparing trees " << std::endl;
  
  //------------------------------------------------------------------------------------------------
  // Read in probes data
  //================================================================================================
  
  //
  // set up input tree
  //
  unsigned int runNum, lumiSec, evtNum;   // event ID
  unsigned int npv;                       // number of primary vertices
  unsigned int pass;                      // whether probe passes requirements
  float        npu;                       // mean number of expected pileup
  float        scale1fb;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  float        met;
  int          njets;
  TLorentzVector *tag=0, *probe=0;        // tag, probe 4-vector
  
  TFile *infile = new TFile(infname.c_str());    assert(infile);
  //TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
  TTree *intree = (TTree*)infile->FindObjectAny("Events"); assert(intree);
  intree->SetBranchAddress("runNum",   &runNum);
  intree->SetBranchAddress("lumiSec",  &lumiSec);
  intree->SetBranchAddress("evtNum",   &evtNum);
  intree->SetBranchAddress("npv",      &npv);
  intree->SetBranchAddress("pass",     &pass);
  intree->SetBranchAddress("npu",      &npu);
  intree->SetBranchAddress("scale1fb", &scale1fb);
  intree->SetBranchAddress("mass",     &mass);
  intree->SetBranchAddress("qtag",     &qtag);
  intree->SetBranchAddress("qprobe",   &qprobe);
  intree->SetBranchAddress("tag",      &tag);
  intree->SetBranchAddress("probe",    &probe);
  intree->SetBranchAddress("met",      &met);
  intree->SetBranchAddress("njets",    &njets);

  //
  // set up output trees
  // Note: each bin is assigned a pass tree and a fail tree
  //
  const unsigned int NBINS_PT     = fPtBinEdgesv.size()-1;
  const unsigned int NBINS_ETA    = fEtaBinEdgesv.size()-1;
  const unsigned int NBINS_PHI    = fPhiBinEdgesv.size()-1;
  const unsigned int NBINS_ETAPT  = NBINS_ETA*NBINS_PT;
  const unsigned int NBINS_ETAPHI = NBINS_ETA*NBINS_PHI;
  const unsigned int NBINS_NPV    = fNPVBinEdgesv.size()-1;
  const unsigned int NBINS_NJETS   = fJetsBinEdgesv.size()-1;
  const unsigned int NBINS_MET    = fMETBinEdgesv.size()-1;
  
  char tname[500];
  float wgt;

  for(unsigned int ibin=0; ibin<NBINS_PT; ibin++) {
    sprintf(tname,"passPt_%i",ibin);
    fPassTreePtv.push_back(new TTree(tname,""));
    fPassTreePtv[ibin]->Branch("m",&mass,"m/F");
    fPassTreePtv[ibin]->Branch("w",&wgt, "w/F");
    fPassTreePtv[ibin]->SetDirectory(0);
    sprintf(tname,"failPt_%i",ibin);
    fFailTreePtv.push_back(new TTree(tname,""));
    fFailTreePtv[ibin]->Branch("m",&mass,"m/F");
    fFailTreePtv[ibin]->Branch("w",&wgt, "w/F");
    fFailTreePtv[ibin]->SetDirectory(0);
  }

  for(unsigned int ibin=0; ibin<NBINS_ETA; ibin++) {
    sprintf(tname,"passEta_%i",ibin);
    fPassTreeEtav.push_back(new TTree(tname,""));
    fPassTreeEtav[ibin]->Branch("m",&mass,"m/F");
    fPassTreeEtav[ibin]->Branch("w",&wgt, "w/F");
    fPassTreeEtav[ibin]->SetDirectory(0);
    sprintf(tname,"failEta_%i",ibin);
    fFailTreeEtav.push_back(new TTree(tname,""));
    fFailTreeEtav[ibin]->Branch("m",&mass,"m/F");
    fFailTreeEtav[ibin]->Branch("w",&wgt, "w/F");
    fFailTreeEtav[ibin]->SetDirectory(0);
  }

  for(unsigned int ibin=0; ibin<NBINS_PHI; ibin++) {
    sprintf(tname,"passPhi_%i",ibin);
    fPassTreePhiv.push_back(new TTree(tname,""));
    fPassTreePhiv[ibin]->Branch("m",&mass,"m/F");
    fPassTreePhiv[ibin]->Branch("w",&wgt, "w/F");
    fPassTreePhiv[ibin]->SetDirectory(0);
    sprintf(tname,"failPhi_%i",ibin);
    fFailTreePhiv.push_back(new TTree(tname,""));
    fFailTreePhiv[ibin]->Branch("m",&mass,"m/F");
    fFailTreePhiv[ibin]->Branch("w",&wgt, "w/F");
    fFailTreePhiv[ibin]->SetDirectory(0);
  }

  for(unsigned int ibin=0; ibin<NBINS_ETAPT; ibin++) {
    sprintf(tname,"passEtaPt_%i",ibin);
    fPassTreeEtaPtv.push_back(new TTree(tname,""));
    fPassTreeEtaPtv[ibin]->Branch("m",&mass,"m/F");
    fPassTreeEtaPtv[ibin]->Branch("w",&wgt, "w/F");
    fPassTreeEtaPtv[ibin]->SetDirectory(0);
    sprintf(tname,"failEtaPt_%i",ibin);
    fFailTreeEtaPtv.push_back(new TTree(tname,""));
    fFailTreeEtaPtv[ibin]->Branch("m",&mass,"m/F");
    fFailTreeEtaPtv[ibin]->Branch("w",&wgt, "w/F");
    fFailTreeEtaPtv[ibin]->SetDirectory(0);
  }

  for(unsigned int ibin=0; ibin<NBINS_ETAPHI; ibin++) {
    sprintf(tname,"passEtaPhi_%i",ibin);
    fPassTreeEtaPhiv.push_back(new TTree(tname,""));
    fPassTreeEtaPhiv[ibin]->Branch("m",&mass,"m/F");
    fPassTreeEtaPhiv[ibin]->Branch("w",&wgt, "w/F");
    fPassTreeEtaPhiv[ibin]->SetDirectory(0);
    sprintf(tname,"failEtaPhi_%i",ibin);
    fFailTreeEtaPhiv.push_back(new TTree(tname,""));
    fFailTreeEtaPhiv[ibin]->Branch("m",&mass,"m/F");
    fFailTreeEtaPhiv[ibin]->Branch("w",&wgt, "w/F");
    fFailTreeEtaPhiv[ibin]->SetDirectory(0);
  }

  for(unsigned int ibin=0; ibin<NBINS_NPV; ibin++) {
    sprintf(tname,"passNPV_%i",ibin);
    fPassTreeNPVv.push_back(new TTree(tname,""));
    fPassTreeNPVv[ibin]->Branch("m",&mass,"m/F");
    fPassTreeNPVv[ibin]->Branch("w",&wgt, "w/F");
    fPassTreeNPVv[ibin]->SetDirectory(0);
    sprintf(tname,"failNPV_%i",ibin);
    fFailTreeNPVv.push_back(new TTree(tname,""));
    fFailTreeNPVv[ibin]->Branch("m",&mass,"m/F");
    fFailTreeNPVv[ibin]->Branch("w",&wgt, "w/F");
    fFailTreeNPVv[ibin]->SetDirectory(0);
  }

  for(unsigned int ibin=0; ibin<NBINS_NJETS; ibin++) {
    sprintf(tname,"passJets_%i",ibin);
    fPassTreeJetsv.push_back(new TTree(tname,""));
    fPassTreeJetsv[ibin]->Branch("m",&mass,"m/F");
    fPassTreeJetsv[ibin]->Branch("w",&wgt, "w/F");
    fPassTreeJetsv[ibin]->SetDirectory(0);
    sprintf(tname,"failJets_%i",ibin);
    fFailTreeJetsv.push_back(new TTree(tname,""));
    fFailTreeJetsv[ibin]->Branch("m",&mass,"m/F");
    fFailTreeJetsv[ibin]->Branch("w",&wgt, "w/F");
    fFailTreeJetsv[ibin]->SetDirectory(0);
  }

  for(unsigned int ibin=0; ibin<NBINS_MET; ibin++) {
    sprintf(tname,"passMET_%i",ibin);
    fPassTreeMETv.push_back(new TTree(tname,""));
    fPassTreeMETv[ibin]->Branch("m",&mass,"m/F");
    fPassTreeMETv[ibin]->Branch("w",&wgt, "w/F");
    fPassTreeMETv[ibin]->SetDirectory(0);
    sprintf(tname,"failMET_%i",ibin);
    fFailTreeMETv.push_back(new TTree(tname,""));
    fFailTreeMETv[ibin]->Branch("m",&mass,"m/F");
    fFailTreeMETv[ibin]->Branch("w",&wgt, "w/F");
    fFailTreeMETv[ibin]->SetDirectory(0);
  }

  std::cout << "   [CEffZFitter] Filling output trees " << std::endl;

  //
  // loop over probes
  //
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    if(qprobe*fCharge < 0) continue;
    if(mass < fFitMassLo) continue;
    if(mass > fFitMassHi) continue;
    if(runNum < fRunNumLo) continue;
    if(runNum > fRunNumHi) continue;

    wgt = scale1fb;
    if(puWeights) {
      wgt *= puWeights->GetBinContent(puWeights->FindBin(npu));
    }
    
    //
    // Find bin indices
    //
    int ipt=-1;
    for(unsigned int ibin=0; ibin<NBINS_PT; ibin++)
      if((probe->Pt() >= fPtBinEdgesv[ibin]) && (probe->Pt() < fPtBinEdgesv[ibin+1]))
        ipt = ibin; 
    if(ipt<0) continue;
    
    int ieta=-1;
    for(unsigned int ibin=0; ibin<NBINS_ETA; ibin++) {
      if(fDoAbsEta) {
        assert(fEtaBinEdgesv[ibin]>=0);
	if((fabs(probe->Eta()) >= fEtaBinEdgesv[ibin]) && (fabs(probe->Eta()) < fEtaBinEdgesv[ibin+1]))
	  ieta = ibin;
      } else {
        if((probe->Eta() >= fEtaBinEdgesv[ibin]) && (probe->Eta() < fEtaBinEdgesv[ibin+1]))
	  ieta = ibin;
      }
    }
    if(ieta<0) continue;
    
    int iphi=-1;
    for(unsigned int ibin=0; ibin<NBINS_PHI; ibin++) {
      if(fDoAbsPhi) {
        assert(fPhiBinEdgesv[ibin]>=0);
	if((fabs(probe->Phi()) >= fPhiBinEdgesv[ibin]) && (fabs(probe->Phi()) < fPhiBinEdgesv[ibin+1]))
	  iphi = ibin;
      } else {
        if((probe->Phi() >= fPhiBinEdgesv[ibin]) && (probe->Phi() < fPhiBinEdgesv[ibin+1]))
	  iphi = ibin;
      }
    }
    if(iphi<0) continue;
    
    int inpv=-1;
    for(unsigned int ibin=0; ibin<NBINS_NPV; ibin++)
      if((npv >= fNPVBinEdgesv[ibin]) && (npv < fNPVBinEdgesv[ibin+1]))
        inpv = ibin; 
    if(inpv<0) continue;

    int ijets=-1;
    for(unsigned int ibin=0; ibin<NBINS_NJETS; ibin++)
      if((njets >= fJetsBinEdgesv[ibin]) && (njets < fJetsBinEdgesv[ibin+1]))
        ijets = ibin; 
    if(ijets<0) continue;
    
    int imet=-1;
    for(unsigned int ibin=0; ibin<NBINS_MET; ibin++)
      if((met >= fMETBinEdgesv[ibin]) && (met < fMETBinEdgesv[ibin+1]))
        imet = ibin; 
    if(imet<0) continue;

    //
    // Fill trees
    //
    if(pass) {
      fPassTreePtv[ipt]->Fill();
      fPassTreeEtav[ieta]->Fill();
      fPassTreePhiv[iphi]->Fill();
      fPassTreeEtaPtv[ipt*NBINS_ETA + ieta]->Fill();
      fPassTreeEtaPhiv[iphi*NBINS_ETA + ieta]->Fill();
      fPassTreeNPVv[inpv]->Fill();
      fPassTreeJetsv[ijets]->Fill();
      fPassTreeMETv[imet]->Fill();
    
    } else {
      fFailTreePtv[ipt]->Fill();
      fFailTreeEtav[ieta]->Fill();
      fFailTreePhiv[iphi]->Fill();
      fFailTreeEtaPtv[ipt*NBINS_ETA + ieta]->Fill();
      fFailTreeEtaPhiv[iphi*NBINS_ETA + ieta]->Fill();
      fFailTreeJetsv[ijets]->Fill();
      fFailTreeMETv[imet]->Fill();
    }
  }
  delete infile;
  infile=0, intree=0;
  
  delete pufile;
  pufile=0, puWeights=0;
  
  fPreparedTrees = true;
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::computeEff()
{
  assert(fIsInitialized);
  assert(fPreparedTrees);
  
  std::cout << "   [CEffZFitter] Computing efficiencies..." << std::endl;

  //------------------------------------------------------------------------------------------------
  // Efficiency calculation
  //================================================================================================
  const unsigned int NBINS_PT     = fPtBinEdgesv.size()-1;
  const unsigned int NBINS_ETA    = fEtaBinEdgesv.size()-1;
  const unsigned int NBINS_PHI    = fPhiBinEdgesv.size()-1;
  const unsigned int NBINS_ETAPT  = NBINS_ETA*NBINS_PT;
  const unsigned int NBINS_ETAPHI = NBINS_ETA*NBINS_PHI;
  const unsigned int NBINS_NPV    = fNPVBinEdgesv.size()-1;
  const unsigned int NBINS_NJETS  = fJetsBinEdgesv.size()-1;
  const unsigned int NBINS_MET    = fMETBinEdgesv.size()-1;
  double ptBinEdges[fPtBinEdgesv.size()];   for(unsigned int i=0; i<fPtBinEdgesv.size();  i++) { ptBinEdges[i]  = fPtBinEdgesv[i];  }
  double etaBinEdges[fEtaBinEdgesv.size()]; for(unsigned int i=0; i<fEtaBinEdgesv.size(); i++) { etaBinEdges[i] = fEtaBinEdgesv[i]; }
  double phiBinEdges[fPhiBinEdgesv.size()]; for(unsigned int i=0; i<fPhiBinEdgesv.size(); i++) { phiBinEdges[i] = fPhiBinEdgesv[i]; }
  
  TGraphAsymmErrors *grEffPt=0;
  TGraphAsymmErrors *grEffEta=0;
  TGraphAsymmErrors *grEffPhi=0;
  TGraphAsymmErrors *grEffNPV=0;
  TGraphAsymmErrors *grEffJets=0;
  TGraphAsymmErrors *grEffMET=0;
  TH2D *hEffEtaPt   = new TH2D("hEffEtaPt","",NBINS_ETA,etaBinEdges,NBINS_PT,ptBinEdges);
  TH2D *hErrlEtaPt  = (TH2D*)hEffEtaPt->Clone("hErrlEtaPt");
  TH2D *hErrhEtaPt  = (TH2D*)hEffEtaPt->Clone("hErrhEtaPt");
  TH2D *hEffEtaPhi  = new TH2D("hEffEtaPhi","",NBINS_ETA,etaBinEdges,NBINS_PHI,phiBinEdges);
  TH2D *hErrlEtaPhi = (TH2D*)hEffEtaPhi->Clone("hErrlEtaPhi");
  TH2D *hErrhEtaPhi = (TH2D*)hEffEtaPhi->Clone("hErrhEtaPhi");

  if(fDoPt)     { grEffPt = makeEffGraph(fPtBinEdgesv, fPassTreePtv, fFailTreePtv, "pt");      grEffPt->SetName("grEffPt"); }  
  if(fDoEta)    { grEffEta = makeEffGraph(fEtaBinEdgesv, fPassTreeEtav, fFailTreeEtav, "eta"); grEffEta->SetName("grEffEta"); }
  if(fDoPhi)    { grEffPhi = makeEffGraph(fPhiBinEdgesv, fPassTreePhiv, fFailTreePhiv, "phi"); grEffPhi->SetName("grEffPhi"); }
  if(fDoNPV)    { grEffNPV = makeEffGraph(fNPVBinEdgesv, fPassTreeNPVv, fFailTreeNPVv, "npv"); grEffNPV->SetName("grEffNPV"); }
  if(fDoJets)    { grEffJets = makeEffGraph(fJetsBinEdgesv, fPassTreeJetsv, fFailTreeJetsv, "jets"); grEffJets->SetName("grEffJets"); }
  if(fDoMET)    { grEffMET = makeEffGraph(fMETBinEdgesv, fPassTreeMETv, fFailTreeMETv, "met"); grEffMET->SetName("grEffMET"); }
  if(fDoEtaPt)  { makeEffHist2D(hEffEtaPt, hErrlEtaPt, hErrhEtaPt, fPassTreeEtaPtv, fFailTreeEtaPtv, "etapt"); }
  if(fDoEtaPhi) { makeEffHist2D(hEffEtaPhi, hErrlEtaPhi, hErrhEtaPhi, fPassTreeEtaPhiv, fFailTreeEtaPhiv, "etaphi"); }
  
  //------------------------------------------------------------------------------------------------
  // Output
  //================================================================================================
  
  //
  // output ROOT file
  //
  std::string outfname = fOutputDir + std::string("/eff.root");
  TFile *outfile = new TFile(outfname.c_str(), "RECREATE");
  if(grEffPt)   grEffPt->Write();
  if(grEffEta)  grEffEta->Write();
  if(grEffPhi)  grEffPhi->Write();
  if(grEffNPV)  grEffNPV->Write();
  if(grEffJets) grEffJets->Write();
  if(grEffMET)  grEffMET->Write();
  hEffEtaPt->Write();
  hErrlEtaPt->Write();
  hErrhEtaPt->Write();
  hEffEtaPhi->Write();
  hErrlEtaPhi->Write();
  hErrhEtaPhi->Write();
  outfile->Close();
  
  //
  // text file of summary tables and summary plots
  //
  ofstream txtfile;
  char txtfname[1000];    
  sprintf(txtfname,"%s/summary.txt",fOutputDir.c_str());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
 
  CEffUser1D effpt;
  CEffUser1D effeta;
  CEffUser1D effphi;
  CEffUser1D effnpv;

  CEffUser2D effetapt;
  CEffUser2D effetaphi;

  TCanvas *c = MakeCanvas("c","c",800,600);

  if(grEffPt) {
    effpt.loadEff(grEffPt);
    effpt.printEff(txtfile);
    txtfile << endl;
    effpt.printErrLow(txtfile);
    txtfile << endl;
    effpt.printErrHigh(txtfile);
    txtfile << endl;
    txtfile << endl;

    CPlot plotEffPt("effpt","","probe p_{T} [GeV/c]","#varepsilon");
    plotEffPt.AddGraph(grEffPt,"",kBlack);
    //plotEffPt.SetYRange(0.6, 1.03);
    plotEffPt.SetYRange(-0.03, 1.03);
    plotEffPt.SetXRange(0.9*(fPtBinEdgesv[0]),1.1*(fPtBinEdgesv[NBINS_PT]));
    plotEffPt.SetLogx();
    plotEffPt.Draw(c,true,"png"); 
    plotEffPt.Draw(c,true,"pdf"); 
  }
  
  if(grEffEta) {
    effeta.loadEff(grEffEta);
    effeta.printEff(txtfile);
    txtfile << endl;
    effeta.printErrLow(txtfile);
    txtfile << endl;
    effeta.printErrHigh(txtfile);
    txtfile << endl;
    txtfile << endl;

    CPlot plotEffEta("effeta","","probe #eta","#varepsilon");
    if(fDoAbsEta) plotEffEta.SetXTitle("probe |#eta|");
    plotEffEta.AddGraph(grEffEta,"",kBlack);
    plotEffEta.SetYRange(0.2, 1.04);
    plotEffEta.Draw(c,true,"png");
    plotEffEta.Draw(c,true,"pdf");

    CPlot plotEffEta2("effeta2","","probe #eta","#varepsilon");
    if(fDoAbsEta) plotEffEta2.SetXTitle("probe |#eta|");
    plotEffEta2.AddGraph(grEffEta,"",kBlack);
    plotEffEta2.SetYRange(0.6, 1.03);
    plotEffEta2.Draw(c,true,"png");
    plotEffEta2.Draw(c,true,"pdf");
  }
  
  if(grEffPhi) {
    effphi.loadEff(grEffPhi);
    effphi.printEff(txtfile);
    txtfile << endl;
    effphi.printErrLow(txtfile);
    txtfile << endl;
    effphi.printErrHigh(txtfile);
    txtfile << endl;
    txtfile << endl;

    CPlot plotEffPhi("effphi","","probe #phi","#varepsilon");
    plotEffPhi.AddGraph(grEffPhi,"",kBlack);
    plotEffPhi.SetYRange(0.6, 1.03);
    plotEffPhi.Draw(c,true,"png"); 
    plotEffPhi.Draw(c,true,"pdf"); 
  }
  
  if(grEffNPV) {
    effnpv.loadEff(grEffNPV);
    effnpv.printEff(txtfile);
    txtfile << endl;
    effnpv.printErrLow(txtfile);
    txtfile << endl;
    effnpv.printErrHigh(txtfile);

    CPlot plotEffNPV("effnpv","","N_{PV}","#varepsilon");
    plotEffNPV.AddGraph(grEffNPV,"",kBlack);
    plotEffNPV.SetYRange(0.6, 1.03);
    plotEffNPV.SetXRange(0.9*(fNPVBinEdgesv[0]),1.1*(fNPVBinEdgesv[NBINS_NPV]));
    plotEffNPV.Draw(c,true,"png"); 
    plotEffNPV.Draw(c,true,"pdf"); 
  }

  if(grEffJets) {
    effnpv.loadEff(grEffJets);
    effnpv.printEff(txtfile);
    txtfile << endl;
    effnpv.printErrLow(txtfile);
    txtfile << endl;
    effnpv.printErrHigh(txtfile);

    CPlot plotEffJets("effjets","","n_{jets}","#varepsilon");
    plotEffJets.AddGraph(grEffJets,"",kBlack);
    plotEffJets.SetYRange(0.6, 1.03);
    plotEffJets.SetXRange(0.9*(fJetsBinEdgesv[0]),1.1*(fJetsBinEdgesv[NBINS_NJETS]));
    plotEffJets.Draw(c,true,"png"); 
    plotEffJets.Draw(c,true,"pdf"); 
  }

  if(grEffMET) {
    effnpv.loadEff(grEffMET);
    effnpv.printEff(txtfile);
    txtfile << endl;
    effnpv.printErrLow(txtfile);
    txtfile << endl;
    effnpv.printErrHigh(txtfile);

    CPlot plotEffMET("effnpv","","N_{PV}","#varepsilon");
    plotEffMET.AddGraph(grEffMET,"",kBlack);
    plotEffMET.SetYRange(0.6, 1.03);
    plotEffMET.SetXRange(0.9*(fMETBinEdgesv[0]),1.1*(fMETBinEdgesv[NBINS_MET]));
    plotEffMET.Draw(c,true,"png"); 
    plotEffMET.Draw(c,true,"pdf"); 
  }


  gStyle->SetPaintTextFormat("3.2f");
  bool use_mit_palette=true; 
  if(use_mit_palette) {
    static Int_t  colors[100];
    Double_t Red[3]    = { 1, 138./255., 163/255.};
    Double_t Green[3]  = { 1, 139./255., 31/255.};
    Double_t Blue[3]   = { 1, 140./255., 52/255.};
    Double_t Length[3] = { 0.00, 0.35, 1.00 };
    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,100);
    for (int i=0; i<100; i++) colors[i] = FI+i;
    gStyle->SetPalette(100,colors);
  } else {
    gStyle->SetPalette(1);
  }
  
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.15);

  if(hEffEtaPt->GetEntries()>0) {
    effetapt.loadEff(hEffEtaPt,hErrlEtaPt,hErrhEtaPt);
    effetapt.printEff(txtfile);     txtfile << endl;
    effetapt.printErrLow(txtfile);  txtfile << endl;
    effetapt.printErrHigh(txtfile); txtfile << endl;
    txtfile << endl;

    hEffEtaPt->SetTitleOffset(1.2,"Y");
    if(NBINS_PT>2) { hEffEtaPt->GetYaxis()->SetRangeUser(fPtBinEdgesv[0],fPtBinEdgesv[NBINS_PT]); }
    CPlot plotEffEtaPt("effetapt","","probe #eta","probe p_{T} [GeV/c]");
    if(fDoAbsEta) { plotEffEtaPt.SetXTitle("probe |#eta|"); }
    plotEffEtaPt.AddHist2D(hEffEtaPt,"TEXT COLZ");
    plotEffEtaPt.SetLogy();
    plotEffEtaPt.Draw(c,true,"png");
    plotEffEtaPt.Draw(c,true,"pdf");

    hErrlEtaPt->SetTitleOffset(1.2,"Y");
    if(NBINS_PT>2) { hErrlEtaPt->GetYaxis()->SetRangeUser(fPtBinEdgesv[0],fPtBinEdgesv[NBINS_PT]); }
    CPlot plotErrlEtaPt("errletapt","","probe #eta","probe p_{T} [GeV/c]");
    if(fDoAbsEta) { plotErrlEtaPt.SetXTitle("probe |#eta|"); }
    plotErrlEtaPt.AddHist2D(hErrlEtaPt,"TEXT COLZ");
    plotErrlEtaPt.SetLogy();
    plotErrlEtaPt.Draw(c,true,"png");
    plotErrlEtaPt.Draw(c,true,"pdf");

    hErrhEtaPt->SetTitleOffset(1.2,"Y");
    if(NBINS_PT>2) { hErrhEtaPt->GetYaxis()->SetRangeUser(fPtBinEdgesv[0],fPtBinEdgesv[NBINS_PT]); }
    CPlot plotErrhEtaPt("errhetapt","","probe #eta","probe p_{T} [GeV/c]");
    if(fDoAbsEta) { plotErrhEtaPt.SetXTitle("probe |#eta|"); }
    plotErrhEtaPt.AddHist2D(hErrhEtaPt,"TEXT COLZ");
    plotErrhEtaPt.SetLogy();
    plotErrhEtaPt.Draw(c,true,"png");
    plotErrhEtaPt.Draw(c,true,"pdf");
  }

  if(hEffEtaPhi->GetEntries()>0) {
    effetaphi.loadEff(hEffEtaPhi,hErrlEtaPhi,hErrhEtaPhi);
    effetaphi.printEff(txtfile);     txtfile << endl;
    effetaphi.printErrLow(txtfile);  txtfile << endl;
    effetaphi.printErrHigh(txtfile); txtfile << endl;
    txtfile << endl;

    hEffEtaPhi->SetTitleOffset(1.2,"Y");
    CPlot plotEffEtaPhi("effetaphi","","probe #eta","probe #phi");
    if(fDoAbsEta) { plotEffEtaPhi.SetXTitle("probe |#eta|"); }
    plotEffEtaPhi.AddHist2D(hEffEtaPhi,"TEXTE COLZ");
    plotEffEtaPhi.Draw(c,true,"png");
    plotEffEtaPhi.Draw(c,true,"pdf");

    hErrlEtaPhi->SetTitleOffset(1.2,"Y");
    CPlot plotErrlEtaPhi("errletaphi","","probe #eta","probe #phi");
    plotErrlEtaPhi.AddHist2D(hErrlEtaPhi,"TEXT COLZ");
    if(fDoAbsEta) { plotErrlEtaPhi.SetXTitle("probe |#eta|"); }
    plotErrlEtaPhi.Draw(c,true,"png");
    plotErrlEtaPhi.Draw(c,true,"pdf");

    hErrhEtaPhi->SetTitleOffset(1.2,"Y");
    CPlot plotErrhEtaPhi("errhetaphi","","probe #eta","probe #phi");
    if(fDoAbsEta) { plotErrhEtaPhi.SetXTitle("probe |#eta|"); }
    plotErrhEtaPhi.AddHist2D(hErrhEtaPhi,"TEXT COLZ");
    plotErrhEtaPhi.Draw(c,true,"png");
    plotErrhEtaPhi.Draw(c,true,"pdf");
  }

  txtfile.close();

  
  //
  // HTML pages
  //
  makeHTML();
  makeHTML("pt",     NBINS_PT);
  makeHTML("eta",    NBINS_ETA);
  makeHTML("phi",    NBINS_PHI);
  makeHTML("npv",    NBINS_NPV);  
  makeHTML("jets",   NBINS_NJETS);  
  makeHTML("MET",    NBINS_MET);  
  makeHTML("etapt",  NBINS_ETAPT);
  makeHTML("etaphi", NBINS_ETAPHI);
  
  delete grEffPt;
  delete grEffEta;
  delete grEffPhi;
  delete grEffNPV;
  delete grEffJets;
  delete grEffMET;
  grEffPt=0, grEffEta=0, grEffPhi=0, grEffNPV=0, grEffJets=0, grEffMET=0;
  
  delete hEffEtaPt;
  delete hErrlEtaPt;
  delete hErrhEtaPt;
  delete hEffEtaPhi;
  delete hErrlEtaPhi;
  delete hErrhEtaPhi;
  hEffEtaPt=0, hErrlEtaPt=0, hErrhEtaPt=0, hEffEtaPhi=0, hErrlEtaPhi=0, hErrhEtaPhi=0;  

  // delete temporary templates file
  //if(fSigPass==2 || fSigFail==2) {
  //  string outfile_name = fOutputDir + "_binnedTemplates.root";
  //  remove(outfile_name.c_str());
  //} else if(fSigPass==4 || fSigFail==4) {
  //  string outfile_name = fOutputDir + "_unbinnedTemplates.root";
  //  remove(outfile_name.c_str());
  //}

}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::parseConf(const std::string conf)
{
  std::ifstream ifs;
  ifs.open(conf.c_str());
  assert(ifs.is_open());
  std::string line;
  int state=0;
  int opts[8];
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    
    double edge;
    std::stringstream ss(line);
    if(state==0) {
      ss >> opts[0] >> opts[1] >> opts[2] >> opts[3] >> opts[4] >> opts[5] >> opts[6] >> opts[7];
      fDoPt     = (opts[0]==1);
      fDoEta    = (opts[1]==1);
      fDoPhi    = (opts[2]==1);
      fDoNPV    = (opts[3]==1);
      fDoJets   = (opts[4]==1);
      fDoMET    = (opts[5]==1);
      fDoEtaPt  = (opts[6]==1);
      fDoEtaPhi = (opts[7]==1);
    
    } else {
      ss >> edge;
      if     (state==1) { fPtBinEdgesv.push_back(edge);  }
      else if(state==2) { fDoAbsEta = (int(edge)==1); state++; }
      else if(state==3) { fEtaBinEdgesv.push_back(edge); }
      else if(state==4) { fDoAbsPhi = (int(edge)==1); state++; }
      else if(state==5) { fPhiBinEdgesv.push_back(edge); }
      else if(state==6) { fNPVBinEdgesv.push_back(edge); }
      else if(state==7) { fJetsBinEdgesv.push_back(edge); }
      else if(state==8) { fMETBinEdgesv.push_back(edge); }
    }
  }
  ifs.close();
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::makeBinnedTemplates(const std::string temfname, const int charge, TH1D *puWeights)
{
  std::cout << "   [CEffZFitter] Creating binned templates... "; std::cout.flush();

  char hname[500];
  
  const unsigned int NBINS_PT  = fPtBinEdgesv.size()-1;
  const unsigned int NBINS_ETA = fEtaBinEdgesv.size()-1;
  const unsigned int NBINS_PHI = fPhiBinEdgesv.size()-1;
  const unsigned int NBINS_NPV = fNPVBinEdgesv.size()-1;
  
  TH1D* passPt[NBINS_PT];
  TH1D* failPt[NBINS_PT];
  for(unsigned int ibin=0; ibin<NBINS_PT; ibin++) {
    sprintf(hname,"passpt_%i",ibin);
    passPt[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi);
    passPt[ibin]->SetDirectory(0);
    sprintf(hname,"failpt_%i",ibin);
    failPt[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
    failPt[ibin]->SetDirectory(0);
  }
  
  TH1D* passEta[NBINS_ETA];
  TH1D* failEta[NBINS_ETA];
  for(unsigned int ibin=0; ibin<NBINS_ETA; ibin++) {
    sprintf(hname,"passeta_%i",ibin);
    passEta[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi);
    passEta[ibin]->SetDirectory(0);
    sprintf(hname,"faileta_%i",ibin);
    failEta[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
    failEta[ibin]->SetDirectory(0);
  }
  
  TH1D* passPhi[NBINS_PHI];  
  TH1D* failPhi[NBINS_PHI];
  for(unsigned int ibin=0; ibin<NBINS_PHI; ibin++) {
    sprintf(hname,"passphi_%i",ibin);
    passPhi[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi);
    passPhi[ibin]->SetDirectory(0);
    sprintf(hname,"failphi_%i",ibin);
    failPhi[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
    failPhi[ibin]->SetDirectory(0);
  }
    
  TH1D* passEtaPt[NBINS_ETA*NBINS_PT];  
  TH1D* failEtaPt[NBINS_ETA*NBINS_PT];
  for(unsigned int ibin=0; ibin<(NBINS_ETA*NBINS_PT); ibin++) {
    sprintf(hname,"passetapt_%i",ibin);
    passEtaPt[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi);
    passEtaPt[ibin]->SetDirectory(0);
    sprintf(hname,"failetapt_%i",ibin);
    failEtaPt[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
    failEtaPt[ibin]->SetDirectory(0);
  }
  
  TH1D* passEtaPhi[NBINS_ETA*NBINS_PHI];
  TH1D* failEtaPhi[NBINS_ETA*NBINS_PHI];
  for(unsigned int ibin=0; ibin<(NBINS_ETA*NBINS_PHI); ibin++) {
    sprintf(hname,"passetaphi_%i",ibin); 
    passEtaPhi[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi);
    passEtaPhi[ibin]->SetDirectory(0);
    sprintf(hname,"failetaphi_%i",ibin);
    failEtaPhi[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
    failEtaPhi[ibin]->SetDirectory(0);
  }

  TH1D* passNPV[NBINS_NPV];  
  TH1D* failNPV[NBINS_NPV];
  for(unsigned int ibin=0; ibin<NBINS_NPV; ibin++) {
    sprintf(hname,"passnpv_%i",ibin);
    passNPV[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi);
    passNPV[ibin]->SetDirectory(0);
    sprintf(hname,"failnpv_%i",ibin);
    failNPV[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
    failNPV[ibin]->SetDirectory(0);
  }
    
  TH1D* passJets[NBINS_NJETS];  
  TH1D* failJets[NBINS_NJETS];
  for(unsigned int ibin=0; ibin<NBINS_NJETS; ibin++) {
    sprintf(hname,"passjets_%i",ibin);
    passJets[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi);
    passJets[ibin]->SetDirectory(0);
    sprintf(hname,"failjets_%i",ibin);
    failJets[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
    failJets[ibin]->SetDirectory(0);
  }
    
  TH1D* passMET[NBINS_MET];  
  TH1D* failMET[NBINS_MET];
  for(unsigned int ibin=0; ibin<NBINS_MET; ibin++) {
    sprintf(hname,"passmet_%i",ibin);
    passMET[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi);
    passMET[ibin]->SetDirectory(0);
    sprintf(hname,"failmet_%i",ibin);
    failMET[ibin] = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
    failMET[ibin]->SetDirectory(0);
  }
    
  unsigned int runNum, lumiSec, evtNum;   // event ID
  unsigned int npv;                       // number of primary vertices
  unsigned int pass;                      // whether probe passes requirements
  float        npu;                       // mean number of expected pileup
  float        scale1fb;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  float        met;
  int          njets;
  TLorentzVector *tag=0, *probe=0;        // tag, probe 4-vector
  
  TFile *infile = new TFile(infname.c_str());    assert(infile);
  //TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
  TTree *intree = (TTree*)infile->FindObjectAny("Events"); assert(intree);
  intree->SetBranchAddress("runNum",   &runNum);
  intree->SetBranchAddress("lumiSec",  &lumiSec);
  intree->SetBranchAddress("evtNum",   &evtNum);
  intree->SetBranchAddress("npv",      &npv);
  intree->SetBranchAddress("pass",     &pass);
  intree->SetBranchAddress("npu",      &npu);
  intree->SetBranchAddress("scale1fb", &scale1fb);
  intree->SetBranchAddress("mass",     &mass);
  intree->SetBranchAddress("qtag",     &qtag);
  intree->SetBranchAddress("qprobe",   &qprobe);
  intree->SetBranchAddress("tag",      &tag);
  intree->SetBranchAddress("probe",    &probe);
  intree->SetBranchAddress("met",      &met);
  intree->SetBranchAddress("njets",    &njets);
  
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    double puWgt = scale1fb;
    if(puWeights)
      puWgt *= puWeights->GetBinContent(puWeights->FindBin(npu));
    
    if(qprobe*charge < 0) continue;
    
    //
    // Find bin indices
    //
    int ipt=-1;
    for(unsigned int ibin=0; ibin<NBINS_PT; ibin++)
      if((probe->Pt() >= fPtBinEdgesv[ibin]) && (probe->Pt() < fPtBinEdgesv[ibin+1]))
        ipt = ibin; 
    if(ipt<0) continue;
    
    int ieta=-1;
    for(unsigned int ibin=0; ibin<NBINS_ETA; ibin++) {
      if(fDoAbsEta) {
        assert(fEtaBinEdgesv[ibin]>=0);
	if((fabs(probe->Eta()) >= fEtaBinEdgesv[ibin]) && (fabs(probe->Eta()) < fEtaBinEdgesv[ibin+1]))
	  ieta = ibin;
      } else {
        if((probe->Eta() >= fEtaBinEdgesv[ibin]) && (probe->Eta() < fEtaBinEdgesv[ibin+1]))
	  ieta = ibin;
      }
    }
    if(ieta<0) continue;
    
    int iphi=-1;
    for(unsigned int ibin=0; ibin<NBINS_PHI; ibin++) {
      if(fDoAbsPhi) {
        assert(fPhiBinEdgesv[ibin]>=0);
	if((fabs(probe->Phi()) >= fPhiBinEdgesv[ibin]) && (fabs(probe->Phi()) < fPhiBinEdgesv[ibin+1]))
	  iphi = ibin;
      } else {
        if((probe->Phi() >= fPhiBinEdgesv[ibin]) && (probe->Phi() < fPhiBinEdgesv[ibin+1]))
	  iphi = ibin;
      }
    }
    if(iphi<0) continue;
    
    int inpv=-1;
    for(unsigned int ibin=0; ibin<NBINS_NPV; ibin++) {
      if((npv >= fNPVBinEdgesv[ibin]) && (npv < fNPVBinEdgesv[ibin+1]))
        inpv = ibin;
    }
    if(inpv<0) continue;
        
    int ijets=-1;
    for(unsigned int ibin=0; ibin<NBINS_NJETS; ibin++)
      if((njets >= fJetsBinEdgesv[ibin]) && (njets < fJetsBinEdgesv[ibin+1]))
        ijets = ibin; 
    if(ijets<0) continue;
    
    int imet=-1;
    for(unsigned int ibin=0; ibin<NBINS_MET; ibin++)
      if((met >= fMETBinEdgesv[ibin]) && (met < fMETBinEdgesv[ibin+1]))
        imet = ibin; 
    if(imet<0) continue;

    //
    // Fill trees
    //
    if(pass) {
      passPt[ipt]->Fill(mass,puWgt);
      passEta[ieta]->Fill(mass,puWgt);
      passPhi[iphi]->Fill(mass,puWgt);
      passEtaPt[ipt*NBINS_ETA + ieta]->Fill(mass,puWgt);
      passEtaPhi[iphi*NBINS_ETA + ieta]->Fill(mass,puWgt);
      passNPV[inpv]->Fill(mass,puWgt);
      passJets[ijets]->Fill(mass,puWgt);
      passMET[imet]->Fill(mass,puWgt);
    
    } else {
      failPt[ipt]->Fill(mass,puWgt);
      failEta[ieta]->Fill(mass,puWgt);
      failPhi[iphi]->Fill(mass,puWgt);
      failEtaPt[ipt*NBINS_ETA + ieta]->Fill(mass,puWgt);
      failEtaPhi[iphi*NBINS_ETA + ieta]->Fill(mass,puWgt);
      failJets[ijets]->Fill(mass,puWgt);
      failMET[imet]->Fill(mass,puWgt);
    }
        
  }
  infile->Close();
  gSystem->mkdir(("templates/"+fOutputDir).c_str(),true);
  string outfile_name = "templates/" + fOutputDir + "/binnedTemplates.root";
  TFile outfile(outfile_name.c_str(), "RECREATE");
  for(unsigned int ibin=0; ibin<NBINS_PT; ibin++) {
    passPt[ibin]->Write();
    failPt[ibin]->Write();
    delete passPt[ibin];
    delete failPt[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_ETA; ibin++) { 
    passEta[ibin]->Write();
    failEta[ibin]->Write();
    delete passEta[ibin];
    delete failEta[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_PHI; ibin++) {
    passPhi[ibin]->Write();
    failPhi[ibin]->Write();
    delete passPhi[ibin];
    delete failPhi[ibin];
  }
  for(unsigned int ibin=0; ibin<(NBINS_ETA*NBINS_PT); ibin++) {
    passEtaPt[ibin]->Write();
    failEtaPt[ibin]->Write();
    delete passEtaPt[ibin];
    delete failEtaPt[ibin];
  }
  for(unsigned int ibin=0; ibin<(NBINS_ETA*NBINS_PHI); ibin++) {
    passEtaPhi[ibin]->Write();
    failEtaPhi[ibin]->Write();
    delete passEtaPhi[ibin];
    delete failEtaPhi[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_NPV; ibin++) {
    passNPV[ibin]->Write();
    failNPV[ibin]->Write();
    delete passNPV[ibin];
    delete failNPV[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_NJETS; ibin++) {
    passJets[ibin]->Write();
    failJets[ibin]->Write();
    delete passJets[ibin];
    delete failJets[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_MET; ibin++) {
    passMET[ibin]->Write();
    failMET[ibin]->Write();
    delete passMET[ibin];
    delete failMET[ibin];
  }
  outfile.Write();
  outfile.Close(); 

  cout << "Done!" << endl;
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::makeUnbinnedTemplates(const std::string temfname, const int charge)
{
  cout << "   [CEffZFitter] Creating unbinned templates... "; cout.flush();
  
  unsigned int runNum, lumiSec, evtNum;   // event ID
  unsigned int npv;                       // number of primary vertices
  unsigned int pass;                      // whether probe passes requirements
  float        npu;                       // mean number of expected pileup
  float        scale1fb;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  float        met;
  int          njets;
  TLorentzVector *tag=0, *probe=0;        // tag, probe 4-vector
  
  char tname[500];
  
  const unsigned int NBINS_PT    = fPtBinEdgesv.size()-1;
  const unsigned int NBINS_ETA   = fEtaBinEdgesv.size()-1;
  const unsigned int NBINS_PHI   = fPhiBinEdgesv.size()-1;
  const unsigned int NBINS_NPV   = fNPVBinEdgesv.size()-1;
  const unsigned int NBINS_NJETS = fJetsBinEdgesv.size()-1;
  const unsigned int NBINS_MET   = fMETBinEdgesv.size()-1;
  
  TTree* passPt[NBINS_PT];
  TTree* failPt[NBINS_PT];
  for(unsigned int ibin=0; ibin<NBINS_PT; ibin++) {
    sprintf(tname,"passpt_%i",ibin);
    passPt[ibin] = new TTree(tname,"");
    passPt[ibin]->Branch("m",&mass,"m/F");
    passPt[ibin]->SetDirectory(0);    
    sprintf(tname,"failpt_%i",ibin);
    failPt[ibin] = new TTree(tname,"");
    failPt[ibin]->Branch("m",&mass,"m/F");
    failPt[ibin]->SetDirectory(0);
  }
  
  TTree* passEta[NBINS_ETA];
  TTree* failEta[NBINS_ETA];
  for(unsigned int ibin=0; ibin<NBINS_ETA; ibin++) {
    sprintf(tname,"passeta_%i",ibin);
    passEta[ibin] = new TTree(tname,"");
    passEta[ibin]->Branch("m",&mass,"m/F");
    passEta[ibin]->SetDirectory(0); 
    sprintf(tname,"faileta_%i",ibin);
    failEta[ibin] = new TTree(tname,"");
    failEta[ibin]->Branch("m",&mass,"m/F");
    failEta[ibin]->SetDirectory(0); 
  }
  
  TTree* passPhi[NBINS_PHI];  
  TTree* failPhi[NBINS_PHI];
  for(unsigned int ibin=0; ibin<NBINS_PHI; ibin++) {
    sprintf(tname,"passphi_%i",ibin);
    passPhi[ibin] = new TTree(tname,"");
    passPhi[ibin]->Branch("m",&mass,"m/F");
    passPhi[ibin]->SetDirectory(0); 
    sprintf(tname,"failphi_%i",ibin);
    failPhi[ibin] = new TTree(tname,"");
    failPhi[ibin]->Branch("m",&mass,"m/F");
    failPhi[ibin]->SetDirectory(0);
  }
    
  TTree* passEtaPt[NBINS_ETA*NBINS_PT];  
  TTree* failEtaPt[NBINS_ETA*NBINS_PT];
  for(unsigned int ibin=0; ibin<(NBINS_ETA*NBINS_PT); ibin++) {
    sprintf(tname,"passetapt_%i",ibin);
    passEtaPt[ibin] = new TTree(tname,"");
    passEtaPt[ibin]->Branch("m",&mass,"m/F");
    passEtaPt[ibin]->SetDirectory(0); 
    sprintf(tname,"failetapt_%i",ibin);
    failEtaPt[ibin] = new TTree(tname,"");
    failEtaPt[ibin]->Branch("m",&mass,"m/F");
    failEtaPt[ibin]->SetDirectory(0);
  }
  
  TTree* passEtaPhi[NBINS_ETA*NBINS_PHI];
  TTree* failEtaPhi[NBINS_ETA*NBINS_PHI];
  for(unsigned int ibin=0; ibin<(NBINS_ETA*NBINS_PHI); ibin++) {
    sprintf(tname,"passetaphi_%i",ibin); 
    passEtaPhi[ibin] = new TTree(tname,"");
    passEtaPhi[ibin]->Branch("m",&mass,"m/F");
    passEtaPhi[ibin]->SetDirectory(0); 
    sprintf(tname,"failetaphi_%i",ibin);
    failEtaPhi[ibin] = new TTree(tname,"");
    failEtaPhi[ibin]->Branch("m",&mass,"m/F");
    failEtaPhi[ibin]->SetDirectory(0);
  }

  TTree* passNPV[NBINS_NPV];  
  TTree* failNPV[NBINS_NPV];
  for(unsigned int ibin=0; ibin<NBINS_NPV; ibin++) {
    sprintf(tname,"passnpv_%i",ibin);
    passNPV[ibin] = new TTree(tname,"");
    passNPV[ibin]->Branch("m",&mass,"m/F");
    passNPV[ibin]->SetDirectory(0); 
    sprintf(tname,"failnpv_%i",ibin);
    failNPV[ibin] = new TTree(tname,"");
    failNPV[ibin]->Branch("m",&mass,"m/F");
    failNPV[ibin]->SetDirectory(0);
  }    

  TTree* passJets[NBINS_NJETS];  
  TTree* failJets[NBINS_NJETS];
  for(unsigned int ibin=0; ibin<NBINS_Jets; ibin++) {
    sprintf(tname,"passjets_%i",ibin);
    passJets[ibin] = new TTree(tname,"");
    passJets[ibin]->Branch("m",&mass,"m/F");
    passJets[ibin]->SetDirectory(0); 
    sprintf(tname,"failjets_%i",ibin);
    failJets[ibin] = new TTree(tname,"");
    failJets[ibin]->Branch("m",&mass,"m/F");
    failJets[ibin]->SetDirectory(0);
  }    
  
  TTree* passMET[NBINS_MET];  
  TTree* failMET[NBINS_MET];
  for(unsigned int ibin=0; ibin<NBINS_MET; ibin++) {
    sprintf(tname,"passmet_%i",ibin);
    passMET[ibin] = new TTree(tname,"");
    passMET[ibin]->Branch("m",&mass,"m/F");
    passMET[ibin]->SetDirectory(0); 
    sprintf(tname,"failmet_%i",ibin);
    failMET[ibin] = new TTree(tname,"");
    failMET[ibin]->Branch("m",&mass,"m/F");
    failMET[ibin]->SetDirectory(0);
  }    
  
  TFile *infile = new TFile(infname.c_str());    assert(infile);
  //TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
  TTree *intree = (TTree*)infile->FindObjectAny("Events"); assert(intree);
  intree->SetBranchAddress("runNum",   &runNum);
  intree->SetBranchAddress("lumiSec",  &lumiSec);
  intree->SetBranchAddress("evtNum",   &evtNum);
  intree->SetBranchAddress("npv",      &npv);
  intree->SetBranchAddress("pass",     &pass);
  intree->SetBranchAddress("npu",      &npu);
  intree->SetBranchAddress("scale1fb", &scale1fb);
  intree->SetBranchAddress("mass",     &mass);
  intree->SetBranchAddress("qtag",     &qtag);
  intree->SetBranchAddress("qprobe",   &qprobe);
  intree->SetBranchAddress("tag",      &tag);
  intree->SetBranchAddress("probe",    &probe);
  intree->SetBranchAddress("met",      &met);
  intree->SetBranchAddress("njets",    &njets);
  
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    if(qprobe*charge < 0) continue;
    if(mass < fFitMassLo) continue;
    if(mass > fFitMassHi) continue;
    
    //
    // Find bin indices
    //
    int ipt=-1;
    for(unsigned int ibin=0; ibin<NBINS_PT; ibin++)
      if((probe->Pt() >= fPtBinEdgesv[ibin]) && (probe->Pt() < fPtBinEdgesv[ibin+1]))
        ipt = ibin; 
    if(ipt<0) continue;
    
    int ieta=-1;
    for(unsigned int ibin=0; ibin<NBINS_ETA; ibin++) {
      if(fDoAbsEta) {
        assert(fEtaBinEdgesv[ibin]>=0);
	if((fabs(probe->Eta()) >= fEtaBinEdgesv[ibin]) && (fabs(probe->Eta()) < fEtaBinEdgesv[ibin+1]))
	  ieta = ibin;
      } else {
        if((probe->Eta() >= fEtaBinEdgesv[ibin]) && (probe->Eta() < fEtaBinEdgesv[ibin+1]))
	  ieta = ibin;
      }
    }
    if(ieta<0) continue;
    
    int iphi=-1;
    for(unsigned int ibin=0; ibin<NBINS_PHI; ibin++) {
      if(fDoAbsPhi) {
        assert(fPhiBinEdgesv[ibin]>=0);
	if((fabs(probe->Phi()) >= fPhiBinEdgesv[ibin]) && (fabs(probe->Phi()) < fPhiBinEdgesv[ibin+1]))
	  iphi = ibin;
      } else {
        if((probe->Phi() >= fPhiBinEdgesv[ibin]) && (probe->Phi() < fPhiBinEdgesv[ibin+1]))
	  iphi = ibin;
      }
    }
    if(iphi<0) continue;
    
    int inpv=-1;
    for(unsigned int ibin=0; ibin<NBINS_NPV; ibin++)
      if((npv >= fNPVBinEdgesv[ibin]) && (npv < fNPVBinEdgesv[ibin+1]))
        inpv = ibin; 
    if(inpv<0) continue;

    int ijets=-1;
    for(unsigned int ibin=0; ibin<NBINS_NJETS; ibin++)
      if((njets >= fJetsBinEdgesv[ibin]) && (njets < fJetsBinEdgesv[ibin+1]))
        ijets = ibin; 
    if(ijets<0) continue;
    
    int imet=-1;
    for(unsigned int ibin=0; ibin<NBINS_MET; ibin++)
      if((met >= fMETBinEdgesv[ibin]) && (met < fMETBinEdgesv[ibin+1]))
        imet = ibin; 
    if(imet<0) continue;
        
    if(pass) {
      passPt[ipt]->Fill();
      passEta[ieta]->Fill();
      passPhi[iphi]->Fill();
      passEtaPt[ipt*NBINS_ETA + ieta]->Fill();
      passEtaPhi[iphi*NBINS_ETA + ieta]->Fill();
      passNPV[inpv]->Fill();
      passJets[ijets]->Fill();
      passMET[imet]->Fill();
    } else {
      failPt[ipt]->Fill();
      failEta[ieta]->Fill();
      failPhi[iphi]->Fill();
      failEtaPt[ipt*NBINS_ETA + ieta]->Fill();
      failEtaPhi[iphi*NBINS_ETA + ieta]->Fill();
      failNPV[inpv]->Fill();
      failJets[ijets]->Fill();
      failMET[imet]->Fill();
    }    
  }
  infile->Close();
  
  string outfile_name = fOutputDir + "_unbinnedTemplates.root";
  TFile outfile(outfile_name.c_str(), "RECREATE");
  for(unsigned int ibin=0; ibin<NBINS_PT; ibin++) {
    passPt[ibin]->Write();
    failPt[ibin]->Write();
    delete passPt[ibin];
    delete failPt[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_ETA; ibin++) { 
    passEta[ibin]->Write();
    failEta[ibin]->Write();
    delete passEta[ibin];
    delete failEta[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_PHI; ibin++) {
    passPhi[ibin]->Write();
    failPhi[ibin]->Write();
    delete passPhi[ibin];
    delete failPhi[ibin];
  }
  for(unsigned int ibin=0; ibin<(NBINS_ETA*NBINS_PT); ibin++) {
    passEtaPt[ibin]->Write();
    failEtaPt[ibin]->Write();
    delete passEtaPt[ibin];
    delete failEtaPt[ibin];
  }
  for(unsigned int ibin=0; ibin<(NBINS_ETA*NBINS_PHI); ibin++) {
    passEtaPhi[ibin]->Write();
    failEtaPhi[ibin]->Write();
    delete passEtaPhi[ibin];
    delete failEtaPhi[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_NPV; ibin++) {
    passNPV[ibin]->Write();
    failNPV[ibin]->Write();
    delete passNPV[ibin];
    delete failNPV[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_NJETS; ibin++) {
    passJets[ibin]->Write();
    failJets[ibin]->Write();
    delete passJets[ibin];
    delete failJets[ibin];
  }
  for(unsigned int ibin=0; ibin<NBINS_MET; ibin++) {
    passMET[ibin]->Write();
    failMET[ibin]->Write();
    delete passMET[ibin];
    delete failMET[ibin];
  }
  outfile.Write();
  outfile.Close(); 

  cout << "Done!" << endl;
}

//--------------------------------------------------------------------------------------------------
TGraphAsymmErrors* CEffZFitter::makeEffGraph(const std::vector<double> &edgesv, 
                                             const std::vector<TTree*> &passv,
					     const std::vector<TTree*> &failv,
                                             const std::string name)
{
  const unsigned n = edgesv.size()-1;
  double xval[n], xerr[n];
  double yval[n], yerrl[n], yerrh[n];
  
  TCanvas *cpass = MakeCanvas("cpass","cpass",720,540);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  cpass->SetTickx(1);
  cpass->SetTicky(1);
  TCanvas *cfail = MakeCanvas("cfail","cfail",720,540); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,cpass->GetWindowTopX()+cfail->GetBorderSize()+540);
  cfail->SetTickx(1);
  cfail->SetTicky(1);  
  
  for(unsigned int ibin=0; ibin<n; ibin++) {
    xval[ibin] = 0.5*(edgesv[ibin+1] + edgesv[ibin]);
    xerr[ibin] = 0.5*(edgesv[ibin+1] - edgesv[ibin]);
    
    double eff, errl, errh;
    
    if(fSigPass==CSignalModel::kNone && fSigFail==CSignalModel::kNone) {  // Cut-and-count
      performCount(eff, errl, errh,
                   ibin,
		   edgesv[ibin], edgesv[ibin+1],
		   0, 0,
                   passv[ibin], failv[ibin], 
                   name, cpass, cfail);
		 
    } else {  // Fit Z peak
      ifstream rfile;
      char rname[512];
      sprintf(rname,"%s/plots/fitres%s_%i.txt",fOutputDir.c_str(),name.c_str(),ibin);
      rfile.open(rname);

      if(rfile.is_open()) {       
        parseFitResults(rfile,eff,errl,errh);
        rfile.close();
        
      } else {
        performFit(eff, errl, errh,
	           ibin, 
		   edgesv[ibin], edgesv[ibin+1], 
		   0, 0,
                   passv[ibin], failv[ibin], 
                   name, cpass, cfail, 4); 
      }
    }  
    
    yval[ibin]  = eff;
    yerrl[ibin] = errl;
    yerrh[ibin] = errh;
  }
  delete cpass;
  delete cfail;
  
  return new TGraphAsymmErrors(n,xval,yval,xerr,xerr,yerrl,yerrh);
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::makeEffHist2D(TH2D *hEff, TH2D *hErrl, TH2D *hErrh, 
                                const std::vector<TTree*> &passv,
		                const std::vector<TTree*> &failv,
                                const std::string name)
{
  TCanvas *cpass = MakeCanvas("cpass","cpass",720,540);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  cpass->SetTickx(1);
  cpass->SetTicky(1);
  TCanvas *cfail = MakeCanvas("cfail","cfail",720,540); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,cpass->GetWindowTopX()+cfail->GetBorderSize()+540);
  cfail->SetTickx(1);
  cfail->SetTicky(1);
  
  for(int iy=0; iy<hEff->GetNbinsY(); iy++) {
    for(int ix=0; ix<hEff->GetNbinsX(); ix++) {
      int ibin = iy*(hEff->GetNbinsX()) + ix;
      
      double eff, errl, errh;

      if(fSigPass==CSignalModel::kNone && fSigFail==CSignalModel::kNone) {  // Cut-and-count
        performCount(eff, errl, errh, 
                     ibin, 
                     hEff->GetXaxis()->GetBinLowEdge(ix+1), hEff->GetXaxis()->GetBinLowEdge(ix+2),
                     hEff->GetYaxis()->GetBinLowEdge(iy+1), hEff->GetYaxis()->GetBinLowEdge(iy+2),
                     passv[ibin], failv[ibin],
                     name, cpass, cfail);
      
      } else {  // Fit Z peak
        ifstream rfile;
        char rname[512];
        sprintf(rname,"%s/plots/fitres%s_%i.txt",fOutputDir.c_str(),name.c_str(),ibin);
        rfile.open(rname);

        if(rfile.is_open()) {
          parseFitResults(rfile,eff,errl,errh);
          rfile.close();

        } else {
          performFit(eff, errl, errh,
                     ibin,
                     hEff->GetXaxis()->GetBinLowEdge(ix+1), hEff->GetXaxis()->GetBinLowEdge(ix+2),
                     hEff->GetYaxis()->GetBinLowEdge(iy+1), hEff->GetYaxis()->GetBinLowEdge(iy+2),
                     passv[ibin], failv[ibin],
                     name, cpass, cfail, 4);
        }
      }

      hEff ->SetBinContent(hEff->GetBin(ix+1, iy+1), eff);
      hEff ->SetBinError(hEff->GetBin(ix+1, iy+1), TMath::Max(errl,errh));
      hErrl->SetBinContent(hErrl->GetBin(ix+1, iy+1), errl);
      hErrh->SetBinContent(hErrh->GetBin(ix+1, iy+1), errh);
    }    
  }
  hEff->SetMaximum(1.);
  hEff->SetMinimum(0.);
  hErrl->SetMaximum(.4);
  hErrl->SetMinimum(0);
  hErrh->SetMaximum(.4);
  hErrh->SetMinimum(0);
  delete cpass;
  delete cfail;  
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::performCount(double &resEff, double &resErrl, double &resErrh,
                               const int ibin,
			       const double xbinLo, const double xbinHi,
			       const double ybinLo, const double ybinHi,
                               TTree *passTree, TTree *failTree, 
                               const std::string name, TCanvas *cpass, TCanvas *cfail)
{
  float m,w;
  char pname[500];
  char binlabelx[1000];
  char binlabely[1000];
  char yield[500];
  char ylabel[500];
  char effstr[1000];    
  
  double npass=0, ntotal=0;
  passTree->SetBranchAddress("m",&m);
  passTree->SetBranchAddress("w",&w);
  for(unsigned int ientry=0; ientry<passTree->GetEntries(); ientry++) {
    passTree->GetEntry(ientry);
    if(m<fMassLo || m>fMassHi) continue;
    npass+=w;
    ntotal+=w;
  }
  failTree->SetBranchAddress("m",&m);
  failTree->SetBranchAddress("w",&w);
  for(unsigned int ientry=0; ientry<failTree->GetEntries(); ientry++) {
    failTree->GetEntry(ientry);
    if(m<fMassLo || m>fMassHi) continue;
    ntotal+=w;
  }
  resEff  = (ntotal>0) ? npass/ntotal : 0;
  if(fUncMethod==0) {
    resErrl = resEff - TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.68269, false);
    resErrh = TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.68269, true) - resEff;
  }
   
  if(name.compare("pt")==0) {
    sprintf(binlabelx,"%i GeV < p_{T} < %i GeV",int(xbinLo),int(xbinHi));
  
  } else if(name.compare("eta")==0) { 
    if(fDoAbsEta) { sprintf(binlabelx,"%.4f < |#eta| < %.4f",xbinLo,xbinHi); }
    else          { sprintf(binlabelx,"%.4f < #eta < %.4f",  xbinLo,xbinHi); }
  
  } else if(name.compare("phi")==0) { 
    if(fDoAbsPhi) { sprintf(binlabelx,"%.4f < |#phi| < %.4f",xbinLo,xbinHi); }
    else          { sprintf(binlabelx,"%.4f < #phi < %.4f",  xbinLo,xbinHi); } 
  
  } else if(name.compare("etapt")==0) {
    if(fDoAbsEta) sprintf(binlabelx,"%.4f < |#eta| < %.4f",xbinLo,xbinHi);
    else         sprintf(binlabelx,"%.4f < #eta < %.4f",xbinLo,xbinHi);    
    sprintf(binlabely,"%i GeV < p_{T} < %i GeV",int(ybinLo),int(ybinHi));
  
  } else if(name.compare("etaphi")==0) {
    if(fDoAbsEta) { sprintf(binlabelx,"%.4f < |#eta| < %.4f",xbinLo,xbinHi); }
    else          { sprintf(binlabelx,"%.4f < #eta < %.4f",  xbinLo,xbinHi); }                                  
    if(fDoAbsPhi) { sprintf(binlabely,"%.4f < |#phi| < %.4f",ybinLo,ybinHi); }
    else          { sprintf(binlabely,"%.4f < #phi < %.4f",  ybinLo,ybinHi); }
  
  } else if(name.compare("npv")==0) { 
    sprintf(binlabelx,"%i #leq N_{PV} < %i",(int)xbinLo,(int)xbinHi);   
  }
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEff,resErrl,resErrh);

  //
  // Plot passing probes
  //
  TH1D *hpass = new TH1D("hpass","",int((fMassHi-fMassLo)/BIN_SIZE_PASS),fMassLo,fMassHi);
  passTree->SetBranchAddress("m",&m);
  passTree->SetBranchAddress("w",&w);
  hpass->Sumw2();
  for(unsigned int ientry=0; ientry<passTree->GetEntries(); ientry++) {
    passTree->GetEntry(ientry);
    hpass->Fill(m,w);
  }
  sprintf(pname,"pass%s_%i",name.c_str(),ibin);
  sprintf(yield,"%i Events",(unsigned int)npass);
  sprintf(ylabel,"Events / %.1f GeV",(double)BIN_SIZE_PASS);
  CPlot plotPass(pname,"","tag-probe mass [GeV]",ylabel);
  plotPass.AddHist1D(hpass,"E");
  plotPass.AddTextBox("Passing probes",0.70,0.93,0.95,0.99,0,kBlack,62,-1);
  plotPass.AddTextBox(binlabelx,0.21,0.84,0.51,0.89,0,kBlack,42,-1);
  if((name.compare("etapt")==0) || (name.compare("etaphi")==0)) {
    plotPass.AddTextBox(binlabely,0.21,0.79,0.51,0.84,0,kBlack,42,-1);        
    plotPass.AddTextBox(yield,0.21,0.74,0.51,0.79,0,kBlack,42,-1);    
  } else {
    plotPass.AddTextBox(yield,0.21,0.79,0.51,0.84,0,kBlack,42,-1);
  }
  plotPass.AddTextBox(effstr,0.69,0.84,0.94,0.89,0,kBlack,42,-1);
  plotPass.Draw(cpass,true,"png");
  plotPass.Draw(cpass,true,"pdf");
  
  //
  // Plot failing probes
  //
  TH1D *hfail = new TH1D("hfail","",int(fMassHi-fMassLo)/BIN_SIZE_FAIL,fMassLo,fMassHi);
  hfail->Sumw2();
  failTree->SetBranchAddress("m",&m);
  failTree->SetBranchAddress("w",&w);
  for(unsigned int ientry=0; ientry<failTree->GetEntries(); ientry++) {
    failTree->GetEntry(ientry);
    hfail->Fill(m,w);
  }
  sprintf(pname,"fail%s_%i",name.c_str(),ibin);
  sprintf(yield,"%i Events",(unsigned int)(ntotal-npass));
  sprintf(ylabel,"Events / %.1f GeV",(double)BIN_SIZE_FAIL);
  CPlot plotFail(pname,"","tag-probe mass [GeV]",ylabel);
  plotFail.AddHist1D(hfail,"E");
  plotFail.AddTextBox("Failing probes",0.70,0.93,0.95,0.99,0,kBlack,62,-1);
  plotFail.AddTextBox(binlabelx,0.21,0.84,0.51,0.89,0,kBlack,42,-1);
  if((name.compare("etapt")==0) || (name.compare("etaphi")==0)) {
    plotFail.AddTextBox(binlabely,0.21,0.79,0.51,0.84,0,kBlack,42,-1);    
    plotFail.AddTextBox(yield,0.21,0.74,0.51,0.79,0,kBlack,42,-1);    
  } else {
    plotFail.AddTextBox(yield,0.21,0.79,0.51,0.84,0,kBlack,42,-1);
  }
  plotFail.AddTextBox(effstr,0.69,0.84,0.94,0.89,0,kBlack,42,-1);
  plotFail.Draw(cfail,true,"png");
  plotFail.Draw(cfail,true,"pdf");
  
  delete hpass;
  delete hfail;
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::performFit(
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
  unsigned int numThreads
) {
  
  char pname[500];
  char binlabelx[1000];
  char binlabely[1000];
  char yield[500];
  char ylabel[500];
  char effstr[1000];
  char nsigstr[1000];
  char nbkgstr[1000];
  char chi2str[1000];
  
  // Define categories
  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);
  
  RooAbsData *dataPass=0;
  RooAbsData *dataFail=0;
  TH1D histPass("histPass","",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi); 
  TH1D histFail("histFail","",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
  RooAbsData *dataCombined=0;
  
  const bool doBinned = true;//(passTree->GetEntries()>1000 && failTree->GetEntries()>1000);
  
  if(doBinned) {
    passTree->Draw("m>>histPass","w");
    failTree->Draw("m>>histFail","w");
    dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),&histPass);
    dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),&histFail);
    //m.setBins(100);  
   
    dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
                                   RooFit::Index(sample),
                                   RooFit::Import("Pass",*((RooDataHist*)dataPass)),
                                   RooFit::Import("Fail",*((RooDataHist*)dataFail)));  
  
  } else {
    dataPass = new RooDataSet("dataPass","dataPass",passTree,RooArgSet(m));
    dataFail = new RooDataSet("dataFail","dataFail",failTree,RooArgSet(m));
    
    dataCombined = new RooDataSet("dataCombined","dataCombined",RooArgList(m),
                                  RooFit::Index(sample),
                                  RooFit::Import("Pass",*((RooDataSet*)dataPass)),
                                  RooFit::Import("Fail",*((RooDataSet*)dataFail))); 
  }
  double ptMax=0, ptMin=8000;
  if(name.compare("pt")==0) {
    ptMax=xbinHi;
    ptMin=xbinLo;
  } else if(name.compare("etapt")==0) {
    ptMax=ybinHi;
    ptMin=ybinLo;
  } 
  createModels(ptMin, ptMax, name, ibin);
   
  // Define free parameters and initial values
  double NsigMax     = doBinned ? histPass.Integral()+histFail.Integral() : passTree->GetEntries()+failTree->GetEntries();
  double NbkgFailMax = doBinned ? histFail.Integral() : failTree->GetEntries();
  double NbkgPassMax = doBinned ? histPass.Integral() : passTree->GetEntries();
  
  Nsig    .removeRange();
  eff     .removeRange();
  NbkgPass.removeRange();
  NbkgFail.removeRange();
  Nsig    .setRange(0,NsigMax);
  eff     .setRange(0,1.0);
  NbkgPass.setRange(0,NbkgPassMax);
  NbkgFail.setRange(0,NbkgFailMax);  
  Nsig    .setVal(ptMin<=30 ? 0.6*NsigMax : 0.8*NsigMax);
  NbkgPass.setVal(0.01*NbkgPassMax);
  
  // cheap way to estimate the efficiency
  if(doBinned)  eff.setVal( 1./ (1. + histFail.GetBinContent(histFail.GetMaximumBin())/histPass.GetBinContent(histPass.GetMaximumBin()) ) );
  else eff.setVal(.8);
  
  if(ptMin<=20) {
    NbkgFail.setVal(0.5*NbkgFailMax);
  } else if(ptMin<=30) {
    NbkgFail.setVal(0.3*NbkgFailMax);  
  } else {
    NbkgFail.setVal(0.05*NbkgFailMax);  
  }

  if(fBkgPass==CBackgroundModel::kNone) NbkgPass.setVal(0);
  if(fBkgFail==CBackgroundModel::kNone) NbkgFail.setVal(0);
    
  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");  
  totalPdf.addPdf(*modelFail,"Fail");

  RooFitResult *fitResult=0;
  if(numThreads<1) numThreads=1;
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::Extended(),
                             RooFit::Strategy(2),
                             RooFit::Minos(RooArgSet(eff)),
                             RooFit::NumCPU(numThreads),
                             RooFit::Save());
 
  // Refit w/o MINOS if MINOS errors are strange...
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5))
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());
  
  resEff  = eff.getVal();
  resErrl = fabs(eff.getErrorLo());
  resErrh = eff.getErrorHi();
    
  if(name.compare("pt")==0) {
    sprintf(binlabelx,"%i GeV < p_{T} < %i GeV",int(xbinLo),int(xbinHi));
  
  } else if(name.compare("eta")==0) { 
    if(fDoAbsEta) { sprintf(binlabelx,"%.4f < |#eta| < %.4f",xbinLo,xbinHi); }
    else          { sprintf(binlabelx,"%.4f < #eta < %.4f",  xbinLo,xbinHi); }
  
  } else if(name.compare("phi")==0) { 
    if(fDoAbsPhi) { sprintf(binlabelx,"%.4f < |#phi| < %.4f",xbinLo,xbinHi); }
    else          { sprintf(binlabelx,"%.4f < #phi < %.4f",  xbinLo,xbinHi); }
  
  } else if(name.compare("etapt")==0) {
    if(fDoAbsEta) { sprintf(binlabelx,"%.4f < |#eta| < %.4f",xbinLo,xbinHi); }
    else          { sprintf(binlabelx,"%.4f < #eta < %.4f",  xbinLo,xbinHi); } 
    sprintf(binlabely,"%i GeV < p_{T} < %i GeV",int(ybinLo),int(ybinHi));
  
  } else if(name.compare("etaphi")==0) {
    if(fDoAbsEta) { sprintf(binlabelx,"%.4f < |#eta| < %.4f",xbinLo,xbinHi); }
    else          { sprintf(binlabelx,"%.4f < #eta < %.4f",  xbinLo,xbinHi); }                                  
    if(fDoAbsPhi) { sprintf(binlabely,"%.4f < |#phi| < %.4f",ybinLo,ybinHi); }
    else          { sprintf(binlabely,"%.4f < #phi < %.4f",  ybinLo,ybinHi); }
  
  } else if(name.compare("npv")==0) { 
    sprintf(binlabelx,"%i #leq N_{PV} < %i",(int)xbinLo,(int)xbinHi); 
  
  } 
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",eff.getVal(),fabs(eff.getErrorLo()),eff.getErrorHi());

  RooPlot *mframePass = m.frame(Bins(int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS));
  dataPass->plotOn(mframePass,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
  modelPass->plotOn(mframePass);
  if(fBkgPass>0)   
    modelPass->plotOn(mframePass,Components("backgroundPass"),LineStyle(kDashed),LineColor(kRed));

  
  RooPlot *mframeFail = m.frame(Bins(int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL));
  dataFail->plotOn(mframeFail,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  modelFail->plotOn(mframeFail);
  modelFail->plotOn(mframeFail,Components("backgroundFail"),LineStyle(kDashed),LineColor(kRed));
  
  if(fFitResDir=="") {
    //
    // Plot passing probes
    //
    sprintf(pname,"pass%s_%i",name.c_str(),ibin);
    sprintf(yield,"%u Events",(int)passTree->GetEntries());
    sprintf(ylabel,"Events / %.1f GeV",(double)BIN_SIZE_PASS);
    sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigPass.getVal(),NsigPass.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/dof = %.3f",mframePass->chiSquare());
    if(fBkgPass!=CBackgroundModel::kNone)
      sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgPass.getVal(),NbkgPass.getPropagatedError(*fitResult));
    CPlot plotPass(pname,mframePass,"","tag-probe mass [GeV]",ylabel);
    plotPass.AddTextBox("Passing probes",0.70,0.93,0.95,0.99,0,kBlack,62,-1);
    plotPass.AddTextBox(binlabelx,0.21,0.84,0.51,0.89,0,kBlack,42,-1);
    if((name.compare("etapt")==0) || (name.compare("etaphi")==0)) {
      plotPass.AddTextBox(binlabely,0.21,0.79,0.51,0.84,0,kBlack,42,-1);        
      plotPass.AddTextBox(yield,0.21,0.74,0.51,0.79,0,kBlack,42,-1);    
    } else {
      plotPass.AddTextBox(yield,0.21,0.79,0.51,0.84,0,kBlack,42,-1);
    }
    plotPass.AddTextBox(effstr,0.69,0.84,0.94,0.89,0,kBlack,42,-1);
    if(fBkgPass!=CBackgroundModel::kNone) {
      plotPass.AddTextBox(0.69,0.68,0.94,0.83,0,kBlack,42,-1,2,nsigstr,nbkgstr);
      //plotPass.AddTextBox(0.69,0.68,0.94,0.83,0,kBlack,42,-1,3,nsigstr,nbkgstr,chi2str);
    } else {
      plotPass.AddTextBox(0.69,0.73,0.94,0.83,0,kBlack,42,-1,1,nsigstr);
      //plotPass.AddTextBox(0.69,0.73,0.94,0.83,0,kBlack,42,-1,2,nsigstr,chi2str);
    }
    plotPass.Draw(cpass,true,"png");
    plotPass.Draw(cpass,true,"pdf");
   
    //
    // Plot failing probes
    //
    sprintf(pname,"fail%s_%i",name.c_str(),ibin);
    sprintf(yield,"%u Events",(int)failTree->GetEntries());
    sprintf(ylabel,"Events / %.1f GeV",(double)BIN_SIZE_FAIL);
    sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigFail.getVal(),NsigFail.getPropagatedError(*fitResult));
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgFail.getVal(),NbkgFail.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/dof = %.3f",mframePass->chiSquare());
    CPlot plotFail(pname,mframeFail,"","tag-probe mass [GeV]",ylabel);
    plotFail.AddTextBox("Failing probes",0.70,0.93,0.95,0.99,0,kBlack,62,-1);
    plotFail.AddTextBox(binlabelx,0.21,0.84,0.51,0.89,0,kBlack,42,-1);
    if((name.compare("etapt")==0) || (name.compare("etaphi")==0)) {
      plotFail.AddTextBox(binlabely,0.21,0.79,0.51,0.84,0,kBlack,42,-1);    
      plotFail.AddTextBox(yield,0.21,0.74,0.51,0.79,0,kBlack,42,-1);    
    } else {
      plotFail.AddTextBox(yield,0.21,0.79,0.51,0.84,0,kBlack,42,-1);
    }
    plotFail.AddTextBox(effstr,0.69,0.84,0.94,0.89,0,kBlack,42,-1);  
    plotFail.AddTextBox(0.69,0.68,0.94,0.83,0,kBlack,42,-1,2,nsigstr,nbkgstr);
    //plotFail.AddTextBox(0.69,0.68,0.94,0.83,0,kBlack,42,-1,3,nsigstr,nbkgstr,chi2str);
    plotFail.Draw(cfail,true,"png");  
    plotFail.Draw(cfail,true,"pdf");  
  } 
  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[1000];    
  if(fFitResDir=="" || fMethodName=="") sprintf(txtfname,"%s/fitres%s_%i.txt",CPlot::sOutDir.Data(),name.c_str(),ibin);
  else sprintf(txtfname, "%s/%s_fitres%s_%i.txt", fFitResDir.c_str(), fMethodName.c_str(), name.c_str(), ibin);
  printf("\n\nSaved results to %s\n\n", txtfname);
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  fitResult->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResult);
  txtfile.close();
  
  //
  // Clean up
  //
  destroyModels();
  delete dataCombined;
  delete dataPass;
  delete dataFail;
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::printCorrelations(std::ostream& os, RooFitResult* res)
{
  ios_base::fmtflags flags = os.flags();
  const RooArgList *parlist = res->correlation("eff");
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(int i=0; i<parlist->getSize(); i++) {
    for(int j=0; j<parlist->getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::parseFitResults(std::ifstream& ifs, double& eff, double& errl, double& errh)
{
  std::string line;
  while(getline(ifs,line)) {
    size_t found = line.find("eff");
    if(found!=string::npos) {
      found = line.find("+/-");
      if(found!=string::npos) {
        std::string varname, initval, pmstr;
        std::stringstream ss(line);
        ss >> varname >> initval >> eff >> pmstr >> errl;
        errh = errl;
        
      } else {
        std::string varname, initval, errstr;
        std::stringstream ss(line);
        ss >> varname >> initval >> eff >> errstr;
        size_t ipos = errstr.find(",");         
        std::string errlstr = errstr.substr(2,ipos-2);
        std::string errhstr = errstr.substr(ipos+2,errstr.length()-ipos-1);
        errl = atof(errlstr.c_str());
        errh = atof(errhstr.c_str());
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::makeHTML()
{
  ofstream htmlfile;
  char htmlfname[1000];
  sprintf(htmlfname,"%s/plots.html",fOutputDir.c_str());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << " PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;

  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effpt.png\"><img src=\"plots/effpt.png\" alt=\"plots/effpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effeta.png\"><img src=\"plots/effeta.png\" alt=\"plots/effeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effeta2.png\"><img src=\"plots/effeta2.png\" alt=\"plots/effeta2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effphi.png\"><img src=\"plots/effphi.png\" alt=\"plots/effphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"pt.html\">pT bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"eta.html\">&eta; bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"phi.html\">&phi; bins</a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effetapt.png\"><img src=\"plots/effetapt.png\" alt=\"plots/effetapt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errletapt.png\"><img src=\"plots/errletapt.png\" alt=\"plots/errletapt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errhetapt.png\"><img src=\"plots/errhetapt.png\" alt=\"plots/errhetapt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"etapt.html\">&eta;-pT bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effetaphi.png\"><img src=\"plots/effetaphi.png\" alt=\"plots/effetaphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errletaphi.png\"><img src=\"plots/errletaphi.png\" alt=\"plots/errletaphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errhetaphi.png\"><img src=\"plots/errhetaphi.png\" alt=\"plots/errhetaphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"etaphi.html\">&eta;-&phi; bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effnpv.png\"><img src=\"plots/effnpv.png\" alt=\"plots/effnpv.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"npv.html\">N_PV bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
}

void CEffZFitter::makeHTML(const std::string name, const unsigned int nbins)
{
  ofstream htmlfile;
  char htmlfname[1000];
  sprintf(htmlfname,"%s/%s.html",fOutputDir.c_str(),name.c_str());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;

  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;    
  unsigned int i;
  for(i=0; i<nbins; i++) {
    if(i%2==0) htmlfile << "<tr>" << endl;    
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pass" << name << "_" << i << ".png\"><img src=\"plots/pass" << name << "_" << i << ".png\"alt=\"plots/pass" << name << "_" << i << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/fail" << name << "_" << i << ".png\"><img src=\"plots/fail" << name << "_" << i << ".png\"alt=\"plots/fail" << name << "_" << i << ".png\" width=\"100%\"></a></td>" << endl;
    if(i%2) htmlfile << "</tr>" << endl;
  }
  if(i%2) {
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
    
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
}

void CEffZFitter::generateToys(
  const int numtoys
) {
  std::cout << "   [CEffZFitter] Generating toys " << std::endl;
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
  assert(fIsInitialized);
  // set up output directory
  gSystem->mkdir(("toys/"+fOutputDir).c_str(),true);

  fSigRefDir = fOutputDir;
  fBkgRefDir = fOutputDir;
      
  TFile *histfile = 0;
  string outfile_name = "templates/" + fOutputDir + "/binnedTemplates.root";
  histfile = new TFile(outfile_name.c_str());
  assert(histfile);

  const unsigned int NBINS_PT     = fPtBinEdgesv.size()-1;
  const unsigned int NBINS_ETA    = fEtaBinEdgesv.size()-1;
  //const unsigned int NBINS_PHI    = fPhiBinEdgesv.size()-1;
  //const unsigned int NBINS_ETAPT  = NBINS_ETA*NBINS_PT;
  //const unsigned int NBINS_ETAPHI = NBINS_ETA*NBINS_PHI;
  //const unsigned int NBINS_NPV    = fNPVBinEdgesv.size()-1;
  double ptBinEdges[fPtBinEdgesv.size()];   for(unsigned int i=0; i<fPtBinEdgesv.size();  i++) { ptBinEdges[i]  = fPtBinEdgesv[i];  }
  //double etaBinEdges[fEtaBinEdgesv.size()]; for(unsigned int i=0; i<fEtaBinEdgesv.size(); i++) { etaBinEdges[i] = fEtaBinEdgesv[i]; }
  //double phiBinEdges[fPhiBinEdgesv.size()]; for(unsigned int i=0; i<fPhiBinEdgesv.size(); i++) { phiBinEdges[i] = fPhiBinEdgesv[i]; }
  
  assert(fDoEtaPt); // only support eta-pt toy study for now
  assert(fSigPass==CSignalModel::kMCTemplateConvGaussianInit); // only support generating toys for MC gaussian binned signal for now
  assert(fSigFail==CSignalModel::kMCTemplateConvGaussianInit); // only support generating toys for MC gaussian binned signal for now
  assert(fBkgPass==CBackgroundModel::kErfcExpoFixed); // only support generating toys for CMSShape background for now
  assert(fBkgFail==CBackgroundModel::kErfcExpoFixed); // only support generating toys for CMSShape background for now

  // Define signal and background models
  // Need to formally put this in a separate function and probably enumerate the models -- to do later
  // For now most of it is yanked from performFit
  for(unsigned int iy=0; iy < NBINS_PT; iy++) { for(unsigned int ix=0; ix < NBINS_ETA; ix++) {
    printf("hi1 \n");
    int ibin = iy*(NBINS_ETA) + ix; 
    double ptMin=ptBinEdges[iy];
    double ptMax=ptBinEdges[iy+1];
    string name="etapt"; 
    char hname[500];
    
    //signal pdfs
    RooRealVar m("m","mass",fFitMassLo,fFitMassHi);
    m.setBins(10000);

    sprintf(hname,"pass%s_%i",name.c_str(),ibin);
    TH1D *hpass = (TH1D*)histfile->Get(hname);
    assert(hpass);
    sigModPass = new CMCTemplateConvGaussian(m,hpass,true, 1, ibin, name, fOutputDir);
    sprintf(hname,"fail%s_%i",name.c_str(),ibin);
    TH1D *hfail = (TH1D*)histfile->Get(hname);
    assert(hfail);
    sigModFail = new CMCTemplateConvGaussian(m,hfail,false, 1, ibin, name, fOutputDir);
    
    //background pdfs
    bkgModPass = new CErfcExpoFixed(m, true, ibin, name, fOutputDir);
    bkgModFail = new CErfcExpoFixed(m, false, ibin, name, fOutputDir);
    
    std::vector<std::string> paramNames = { "Nsig", "eff", "NbkgPass", "NbkgFail"};
    std::vector<double> params = sigModPass->readSigParams(ibin, name, paramNames, fOutputDir);
    Nsig    .removeRange();
    eff     .removeRange();
    NbkgPass.removeRange();
    NbkgFail.removeRange();
    Nsig.setVal(params[0]);
    eff.setVal(params[1]);
    NbkgPass.setVal(params[2]);
    NbkgFail.setVal(params[3]);
    
    createModels(ptMin, ptMax, name, ibin);
    RooMCStudy *passToyMachine = new RooMCStudy(*modelPass, RooArgSet(m));
    RooMCStudy *failToyMachine = new RooMCStudy(*modelFail, RooArgSet(m));
    
    char toyDir[512], toyFileName[512];
    // make directories for the toys
    for(int iToy=0; iToy<numtoys; iToy++) {
      sprintf(toyDir, "toys/%s/toy%06d", fOutputDir.c_str(), iToy);
      gSystem->mkdir(toyDir, true);

    }
    sprintf(toyFileName, "toys/%s/toy%%06d/pass_etapt_%d.dat", fOutputDir.c_str(), ibin);
    passToyMachine->generate(numtoys, params[0]*params[1]+params[2], kFALSE, toyFileName);
    
    sprintf(toyFileName, "toys/%s/toy%%06d/fail_etapt_%d.dat", fOutputDir.c_str(), ibin);
    failToyMachine->generate(numtoys, params[0]*(1.-params[1])+params[3] , kFALSE, toyFileName);
    
    delete passToyMachine,
    delete failToyMachine,
    delete esignalPass;
    delete ebackgroundPass;
    delete esignalFail;
    delete ebackgroundFail;
    delete modelPass;
    delete modelFail;  
    delete sigModPass;
    delete bkgModPass;
    delete sigModFail;
    delete bkgModFail;        
    //delete datfile;
  }}
  delete histfile;
}

void CEffZFitter::fitToy(
  const std::string methodname,
  const int toynum
) {
  fMethodName=methodname;
  assert(fIsInitialized);
  std::cout << "   [CEffZFitter] Fitting toy... " << std::endl;

  const unsigned int NBINS_PT     = fPtBinEdgesv.size()-1;
  const unsigned int NBINS_ETA    = fEtaBinEdgesv.size()-1;
  //const unsigned int NBINS_PHI    = fPhiBinEdgesv.size()-1;
  //const unsigned int NBINS_ETAPT  = NBINS_ETA*NBINS_PT;
  //const unsigned int NBINS_ETAPHI = NBINS_ETA*NBINS_PHI;
  //const unsigned int NBINS_NPV    = fNPVBinEdgesv.size()-1;
  double ptBinEdges[fPtBinEdgesv.size()];   for(unsigned int i=0; i<fPtBinEdgesv.size();  i++) { ptBinEdges[i]  = fPtBinEdgesv[i];  }
  double etaBinEdges[fEtaBinEdgesv.size()]; for(unsigned int i=0; i<fEtaBinEdgesv.size(); i++) { etaBinEdges[i] = fEtaBinEdgesv[i]; }
  //double phiBinEdges[fPhiBinEdgesv.size()]; for(unsigned int i=0; i<fPhiBinEdgesv.size(); i++) { phiBinEdges[i] = fPhiBinEdgesv[i]; }
  
  assert(fDoEtaPt); // only support eta-pt toy study for now

  // Define signal and background models
  // For now most of it is yanked from performFit
  for(unsigned int iy=0; iy < NBINS_PT; iy++) { for(unsigned int ix=0; ix < NBINS_ETA; ix++) {
    int ibin = iy*(NBINS_ETA) + ix; 
    double ptMin=ptBinEdges[iy];
    double ptMax=ptBinEdges[iy+1];
    double etaMin = etaBinEdges[ix];
    double etaMax = etaBinEdges[ix+1];
    string name="etapt"; 
    
    double eff, errl, errh;

    char toyFileName[512];
    TBranch *weightBranch;
    float w=1;
    
    TTree *passTree = new TTree("passTree", "passTree");
    passTree->SetDirectory(0);
    sprintf(toyFileName, "toys/%s/toy%06d/pass_etapt_%d.dat", fOutputDir.c_str(), toynum, ibin);
    passTree->ReadFile(toyFileName, "m");
    weightBranch = passTree->Branch("w", &w, "w/F");
    for(Long64_t i=0; i<passTree->GetEntries(); i++) {
      weightBranch->Fill();
    }
    TTree *failTree = new TTree("failTree", "failTree");
    failTree->SetDirectory(0);
    sprintf(toyFileName, "toys/%s/toy%06d/fail_etapt_%d.dat", fOutputDir.c_str(), toynum, ibin);
    failTree->ReadFile(toyFileName, "m");
    weightBranch = failTree->Branch("w", &w, "w/F");
    for(Long64_t i=0; i<failTree->GetEntries(); i++) {
      weightBranch->Fill();
    }
    
    //garbage canvases
    TCanvas *cpass = new TCanvas ("cpass", "cpass");
    TCanvas *cfail = new TCanvas ("cfail", "cfail");
    char fitResDirStr[512];
    sprintf(fitResDirStr, "toys/%s/toy%06d", fOutputDir.c_str(), toynum);
    fFitResDir = string(fitResDirStr);
    //printf("fFitResDir=\"%s\"\n", fFitResDir.c_str());
    performFit(eff, errl, errh, ibin, etaMin, etaMax, ptMin, ptMax, passTree, failTree, name, cpass, cfail, 1);
    
    sprintf(toyFileName, "toys/%s/toy%06d/eff_%s_etapt_%d.txt", fOutputDir.c_str(), toynum, methodname.c_str(), ibin);
    ofstream efficiencyFile;
    efficiencyFile.open(toyFileName);
    efficiencyFile << eff << " " << errl << " " << errh << std::endl;
    efficiencyFile.close();

    delete passTree;
    delete failTree;
    delete cpass;
    delete cfail;
  }}

}
void CEffZFitter::destroyModels() {
  if(modelPass) delete modelPass;
  if(modelFail) delete modelFail;
  if(esignalPass)     delete esignalPass;
  if(ebackgroundPass) delete ebackgroundPass;
  if(esignalFail)     delete esignalFail;
  if(ebackgroundFail) delete ebackgroundFail;
  if(sigModPass) delete sigModPass;
  if(bkgModPass) delete bkgModPass;
  if(sigModFail) delete sigModFail;
  if(bkgModFail) delete bkgModFail;
}
void CEffZFitter::createModels(
  double ptMin, double ptMax, const std::string name, const int ibin
) {
  assert(fIsInitialized);
  TFile *histfile = 0, *datfile = 0;
  if(
    fSigPass==CSignalModel::kMCTemplateConvGaussian || 
    fSigPass==CSignalModel::kMCTemplateConvGaussianInit ||
    fSigFail==CSignalModel::kMCTemplateConvGaussian || 
    fSigFail==CSignalModel::kMCTemplateConvGaussianInit
  ) {
    string outfile_name = "templates/" + fOutputDir + "/binnedTemplates.root";
    histfile = new TFile(outfile_name.c_str());
    assert(histfile);
  } else if(
    fSigPass==CSignalModel::kMCDatasetConvGaussian ||
    fSigFail==CSignalModel::kMCDatasetConvGaussian
  ) {
    string outfile_name = "templates/" + fOutputDir + "/unbinnedTemplates.root";
    datfile = new TFile(outfile_name.c_str());
    assert(datfile);
  }
 
  // Define signal and background models
  if(fSigPass == CSignalModel::kBreitWignerConvCrystalBall) {
    sigModPass = new CBreitWignerConvCrystalBall(m,true);
  } else if(fSigPass==CSignalModel::kMCTemplateConvGaussian || fSigPass==CSignalModel::kMCTemplateConvGaussianInit) { 
    char hname[500];
    sprintf(hname,"pass%s_%i",name.c_str(),ibin);
    TH1D *h = (TH1D*)histfile->Get(hname);
    assert(h);
    if(fSigPass==CSignalModel::kMCTemplateConvGaussian)
      sigModPass = new CMCTemplateConvGaussian(m,h,true);
    else if(fSigPass==CSignalModel::kMCTemplateConvGaussianInit)
      sigModPass = new CMCTemplateConvGaussian(m,h,true,1,ibin,name,fSigRefDir);
  } else if(fSigPass==CSignalModel::kVoigtianCBShape) {
    sigModPass = new CVoigtianCBShape(m,true);
  } else if(fSigPass==CSignalModel::kMCDatasetConvGaussian) {
    char tname[500];
    sprintf(tname,"pass%s_%i",name.c_str(),ibin);
    TTree *t = (TTree*)datfile->Get(tname);
    assert(t);
    sigModPass = new CMCDatasetConvGaussian(m,t,true);
  } else if(fSigPass==CSignalModel::kBWCBPlusVoigt) {
    sigModPass = new CBWCBPlusVoigt(m, true, ptMin, ptMax);
  } else if(fSigPass==CSignalModel::kBWCBPlusVoigtBounded) {
    sigModPass = new CBWCBPlusVoigt(m, true, ptMin, ptMax, ibin, name, fSigRefDir);
  }

  //Passing background model
  if(fBkgPass==CBackgroundModel::kExponential) { 
    bkgModPass = new CExponential(m,true);
  } else if(fBkgPass==CBackgroundModel::kErfcExpo) {
    bkgModPass = new CErfcExpo(m,true);
  } else if(fBkgPass==CBackgroundModel::kErfcExpoFixed) {
    bkgModPass = new CErfcExpoFixed(m, true, ibin, name, fBkgRefDir);
  } else if(fBkgPass==CBackgroundModel::kDoubleExp) {
    bkgModPass = new CDoubleExp(m,true);
  } else if(fBkgPass==CBackgroundModel::kLinearExp) {
    bkgModPass = new CLinearExp(m,true);
  } else if(fBkgPass==CBackgroundModel::kQuadraticExp) {
    bkgModPass = new CQuadraticExp(m,true);
  }

  //Failing signal model
  if(fSigFail == CSignalModel::kBreitWignerConvCrystalBall) {
    sigModFail = new CBreitWignerConvCrystalBall(m,false);
  } else if(fSigFail==CSignalModel::kMCTemplateConvGaussian || fSigFail==CSignalModel::kMCTemplateConvGaussianInit) { 
    char hname[500];
    sprintf(hname,"fail%s_%i",name.c_str(),ibin);
    TH1D *h = (TH1D*)histfile->Get(hname);
    assert(h);
    if(fSigFail==CSignalModel::kMCTemplateConvGaussian)
      sigModFail = new CMCTemplateConvGaussian(m,h,false);
    else if(fSigFail==CSignalModel::kMCTemplateConvGaussianInit)
      sigModFail = new CMCTemplateConvGaussian(m,h,false,1,ibin,name,fSigRefDir);
  } else if(fSigFail==CSignalModel::kVoigtianCBShape) {
    sigModFail = new CVoigtianCBShape(m,false);
  } else if(fSigFail==CSignalModel::kMCDatasetConvGaussian) {
    char tname[500];
    sprintf(tname,"fail%s_%i",name.c_str(),ibin);
    TTree *t = (TTree*)datfile->Get(tname);
    assert(t);
    sigModFail = new CMCDatasetConvGaussian(m,t,false);
  } else if(fSigFail==CSignalModel::kBWCBPlusVoigt) {
    sigModFail = new CBWCBPlusVoigt(m, false, ptMin, ptMax);
  } else if(fSigFail==CSignalModel::kBWCBPlusVoigtBounded) {
    sigModFail = new CBWCBPlusVoigt(m, false, ptMin, ptMax, ibin, name, fSigRefDir);
  }

  //Failing background model
  if(fBkgFail==CBackgroundModel::kExponential) { 
    bkgModFail = new CExponential(m,false);
  } else if(fBkgFail==CBackgroundModel::kErfcExpo) {
    bkgModFail = new CErfcExpo(m,false);
  } else if(fBkgFail==CBackgroundModel::kErfcExpoFixed) {
    bkgModFail = new CErfcExpoFixed(m, false, ibin, name, fBkgRefDir);
  } else if(fBkgFail==CBackgroundModel::kDoubleExp) {
    bkgModFail = new CDoubleExp(m,false);
  } else if(fBkgFail==CBackgroundModel::kLinearExp) {
    bkgModFail = new CLinearExp(m,false);
  } else if(fBkgFail==CBackgroundModel::kQuadraticExp) {
    bkgModFail = new CQuadraticExp(m,false);
  }

  if(fMassLo!=fFitMassLo || fMassHi!=fFitMassHi) {
      m.setRange("signalRange",fMassLo,fMassHi);
      
      esignalPass     = new RooExtendPdf("esignalPass","esignalPass",*(sigModPass->model),NsigPass,"signalRange");
      ebackgroundPass = new RooExtendPdf("ebackgroundPass","ebackgroundPass",(fBkgPass>0) ? *(bkgModPass->model) : *(sigModPass->model),NbkgPass,"signalRange");
      modelPass       = new RooAddPdf("modelPass","Model for PASS sample",(fBkgPass>0) ? RooArgList(*esignalPass,*ebackgroundPass) : RooArgList(*esignalPass));    
  
      esignalFail     = new RooExtendPdf("esignalFail","esignalFail",*(sigModFail->model),NsigFail,"signalRange");
      ebackgroundFail = new RooExtendPdf("ebackgroundFail","ebackgroundFail",*(bkgModFail->model),NbkgFail,"signalRange");
      modelFail       = new RooAddPdf("modelFail","Model for FAIL sample", (fBkgFail>0) ? RooArgList(*esignalFail,*ebackgroundFail) : RooArgList(*esignalFail));
    
    } else {
      modelPass = new RooAddPdf("modelPass","Model for PASS sample",
                                (fBkgPass!=CBackgroundModel::kNone) ? RooArgList(*(sigModPass->model),*(bkgModPass->model)) :  RooArgList(*(sigModPass->model)),
                                (fBkgPass!=CBackgroundModel::kNone) ? RooArgList(NsigPass,NbkgPass) : RooArgList(NsigPass));
    
      //modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(sigModFail->model),*(bkgModFail->model)),RooArgList(NsigFail,NbkgFail));
      modelFail = new RooAddPdf("modelFail","Model for FAIL sample",
                                (fBkgFail!=CBackgroundModel::kNone) ? RooArgList(*(sigModFail->model),*(bkgModFail->model)) :  RooArgList(*(sigModFail->model)),
                                (fBkgFail!=CBackgroundModel::kNone) ? RooArgList(NsigFail,NbkgFail) : RooArgList(NsigFail));
    }

}
