#include "CEffZFitter.hh"
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>

#include <cassert>
#include <sstream>
#include <iomanip>

#include "CPlot.hh"
#include "KStyle.hh"
#include "ZSignals.hh"
#include "ZBackgrounds.hh"
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

// bin size constants
#define BIN_SIZE_PASS 2
#define BIN_SIZE_FAIL 2

//--------------------------------------------------------------------------------------------------
CEffZFitter::CEffZFitter():
fIsInitialized(false),
fSigPass      (0),
fBkgPass      (0),
fSigFail      (0),
fBkgFail      (0),
fMassLo       (60),
fMassHi       (120),
fFitMassLo    (60),
fFitMassHi    (120),
fUncMethod    (0),
fOutputDir    ("."),
fDoAbsEta     (false),
fDoAbsPhi     (false),
fDoPt         (false),
fDoEta        (false),
fDoPhi        (false),
fDoEtaPt      (false),
fDoEtaPhi     (false),
fDoNPV        (false)
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
  for(unsigned int i=0; i<fFailTreePtv.size(); i++)     { delete fFailTreePtv[i];     fFailTreePtv[i]=0; }
  for(unsigned int i=0; i<fFailTreeEtav.size(); i++)    { delete fFailTreeEtav[i];    fFailTreeEtav[i]=0; }
  for(unsigned int i=0; i<fFailTreePhiv.size(); i++)    { delete fFailTreePhiv[i];    fFailTreePhiv[i]=0; }
  for(unsigned int i=0; i<fFailTreeEtaPtv.size(); i++)  { delete fFailTreeEtaPtv[i];  fFailTreeEtaPtv[i]=0; }
  for(unsigned int i=0; i<fFailTreeEtaPhiv.size(); i++) { delete fFailTreeEtaPhiv[i]; fFailTreeEtaPhiv[i]=0; }
  for(unsigned int i=0; i<fFailTreeNPVv.size(); i++)    { delete fFailTreeNPVv[i];    fFailTreeNPVv[i]=0; }
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::initialize(
  const std::string conf,
  const int sigpass,
  const int bkgpass,
  const int sigfail,
  const int bkgfail,
  const std::string infname,
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
  const unsigned int runNumLo, const unsigned int runNumHi)
{
  std::cout << "   [CEffZFitter] Initializing... " << std::endl;

  //parse binning configuration file
  parseConf(conf);
  
  fSigPass = sigpass;
  fBkgPass = bkgpass;
  fSigFail = sigfail;
  fBkgFail = bkgfail;

  fMassLo    = massLo;
  fMassHi    = massHi;
  fFitMassLo = fitMassLo;
  fFitMassHi = fitMassHi;

  // set up output directory
  fOutputDir = outdir;
  fSigRefDir = sigRefDir;
  fBkgRefDir = bkgRefDir;
  gSystem->mkdir(fOutputDir.c_str(),true);
  CPlot::sOutDir = TString(outdir.c_str()) + TString("/plots");
      
  TFile *pufile=0;
  TH1D *puWeights=0;
  if(pufname.compare("none")!=0) {
    pufile = new TFile(pufname.c_str());      assert(pufile);
    puWeights = (TH1D*)pufile->Get("pileup"); assert(puWeights); 
  }
  
  // generate templates from MC if necessary
  if(fSigPass==2 || fSigFail==2) {
    makeBinnedTemplates(temfname, charge, puWeights);
  } else if(fSigPass==4 || fSigFail==4) {
    makeUnbinnedTemplates(temfname, charge);
  }
  
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

  //
  // loop over probes
  //
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    if(qprobe*charge < 0) continue;
    if(mass < fFitMassLo) continue;
    if(mass > fFitMassHi) continue;
    if(runNum < runNumLo) continue;
    if(runNum > runNumHi) continue;

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
    
    } else {
      fFailTreePtv[ipt]->Fill();
      fFailTreeEtav[ieta]->Fill();
      fFailTreePhiv[iphi]->Fill();
      fFailTreeEtaPtv[ipt*NBINS_ETA + ieta]->Fill();
      fFailTreeEtaPhiv[iphi*NBINS_ETA + ieta]->Fill();
      fFailTreeNPVv[inpv]->Fill();
    }
  }
  delete infile;
  infile=0, intree=0;
  
  delete pufile;
  pufile=0, puWeights=0;
  
  fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::computeEff()
{
  assert(fIsInitialized);
  
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
  double ptBinEdges[fPtBinEdgesv.size()];   for(unsigned int i=0; i<fPtBinEdgesv.size();  i++) { ptBinEdges[i]  = fPtBinEdgesv[i];  }
  double etaBinEdges[fEtaBinEdgesv.size()]; for(unsigned int i=0; i<fEtaBinEdgesv.size(); i++) { etaBinEdges[i] = fEtaBinEdgesv[i]; }
  double phiBinEdges[fPhiBinEdgesv.size()]; for(unsigned int i=0; i<fPhiBinEdgesv.size(); i++) { phiBinEdges[i] = fPhiBinEdgesv[i]; }
  
  TGraphAsymmErrors *grEffPt=0;
  TGraphAsymmErrors *grEffEta=0;
  TGraphAsymmErrors *grEffPhi=0;
  TGraphAsymmErrors *grEffNPV=0;
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
  if(grEffPt)  grEffPt->Write();
  if(grEffEta) grEffEta->Write();
  if(grEffPhi) grEffPhi->Write();
  if(grEffNPV) grEffNPV->Write();
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
    plotEffPt.SetYRange(0.6, 1.03);
    plotEffPt.SetXRange(0.9*(fPtBinEdgesv[0]),1.1*(fPtBinEdgesv[NBINS_PT-1]));
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
    plotEffNPV.SetXRange(0.9*(fNPVBinEdgesv[0]),1.1*(fNPVBinEdgesv[NBINS_NPV-1]));
    plotEffNPV.Draw(c,true,"png"); 
    plotEffNPV.Draw(c,true,"pdf"); 
  }


  gStyle->SetPalette(1);
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.15);

  if(hEffEtaPt->GetEntries()>0) {
    effetapt.loadEff(hEffEtaPt,hErrlEtaPt,hErrhEtaPt);
    effetapt.printEff(txtfile);     txtfile << endl;
    effetapt.printErrLow(txtfile);  txtfile << endl;
    effetapt.printErrHigh(txtfile); txtfile << endl;
    txtfile << endl;

    hEffEtaPt->SetTitleOffset(1.2,"Y");
    if(NBINS_PT>2) { hEffEtaPt->GetYaxis()->SetRangeUser(fPtBinEdgesv[0],fPtBinEdgesv[NBINS_PT-2]); }
    CPlot plotEffEtaPt("effetapt","","probe #eta","probe p_{T} [GeV/c]");
    if(fDoAbsEta) { plotEffEtaPt.SetXTitle("probe |#eta|"); }
    plotEffEtaPt.AddHist2D(hEffEtaPt,"COLZ");
    plotEffEtaPt.Draw(c,true,"png");
    plotEffEtaPt.Draw(c,true,"pdf");

    hErrlEtaPt->SetTitleOffset(1.2,"Y");
    if(NBINS_PT>2) { hErrlEtaPt->GetYaxis()->SetRangeUser(fPtBinEdgesv[0],fPtBinEdgesv[NBINS_PT-2]); }
    CPlot plotErrlEtaPt("errletapt","","probe #eta","probe p_{T} [GeV/c]");
    if(fDoAbsEta) { plotErrlEtaPt.SetXTitle("probe |#eta|"); }
    plotErrlEtaPt.AddHist2D(hErrlEtaPt,"COLZ");
    plotErrlEtaPt.Draw(c,true,"png");
    plotErrlEtaPt.Draw(c,true,"pdf");

    hErrhEtaPt->SetTitleOffset(1.2,"Y");
    if(NBINS_PT>2) { hErrhEtaPt->GetYaxis()->SetRangeUser(fPtBinEdgesv[0],fPtBinEdgesv[NBINS_PT-2]); }
    CPlot plotErrhEtaPt("errhetapt","","probe #eta","probe p_{T} [GeV/c]");
    if(fDoAbsEta) { plotErrhEtaPt.SetXTitle("probe |#eta|"); }
    plotErrhEtaPt.AddHist2D(hErrhEtaPt,"COLZ");
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
    plotEffEtaPhi.AddHist2D(hEffEtaPhi,"COLZ");
    plotEffEtaPhi.Draw(c,true,"png");
    plotEffEtaPhi.Draw(c,true,"pdf");

    hErrlEtaPhi->SetTitleOffset(1.2,"Y");
    CPlot plotErrlEtaPhi("errletaphi","","probe #eta","probe #phi");
    plotErrlEtaPhi.AddHist2D(hErrlEtaPhi,"COLZ");
    if(fDoAbsEta) { plotErrlEtaPhi.SetXTitle("probe |#eta|"); }
    plotErrlEtaPhi.Draw(c,true,"png");
    plotErrlEtaPhi.Draw(c,true,"pdf");

    hErrhEtaPhi->SetTitleOffset(1.2,"Y");
    CPlot plotErrhEtaPhi("errhetaphi","","probe #eta","probe #phi");
    if(fDoAbsEta) { plotErrhEtaPhi.SetXTitle("probe |#eta|"); }
    plotErrhEtaPhi.AddHist2D(hErrhEtaPhi,"COLZ");
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
  makeHTML("etapt",  NBINS_ETAPT);
  makeHTML("etaphi", NBINS_ETAPHI);
  
  delete grEffPt;
  delete grEffEta;
  delete grEffPhi;
  delete grEffNPV;
  grEffPt=0, grEffEta=0, grEffPhi=0, grEffNPV=0;
  
  delete hEffEtaPt;
  delete hErrlEtaPt;
  delete hErrhEtaPt;
  delete hEffEtaPhi;
  delete hErrlEtaPhi;
  delete hErrhEtaPhi;
  hEffEtaPt=0, hErrlEtaPt=0, hErrhEtaPt=0, hEffEtaPhi=0, hErrlEtaPhi=0, hErrhEtaPhi=0;  

  // delete temporary templates file
  if(fSigPass==2 || fSigFail==2) {
    string outfile_name = fOutputDir + "_binnedTemplates.root";
    remove(outfile_name.c_str());
  } else if(fSigPass==4 || fSigFail==4) {
    string outfile_name = fOutputDir + "_unbinnedTemplates.root";
    remove(outfile_name.c_str());
  }

}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::parseConf(const std::string conf)
{
  std::ifstream ifs;
  ifs.open(conf.c_str());
  assert(ifs.is_open());
  std::string line;
  int state=0;
  int opts[6];
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    
    double edge;
    std::stringstream ss(line);
    if(state==0) {
      ss >> opts[0] >> opts[1] >> opts[2] >> opts[3] >> opts[4] >> opts[5];
      fDoPt     = (opts[0]==1);
      fDoEta    = (opts[1]==1);
      fDoPhi    = (opts[2]==1);
      fDoNPV    = (opts[3]==1);
      fDoEtaPt  = (opts[4]==1);
      fDoEtaPhi = (opts[5]==1);
    
    } else {
      ss >> edge;
      if     (state==1) { fPtBinEdgesv.push_back(edge);  }
      else if(state==2) { fDoAbsEta = (int(edge)==1); state++; }
      else if(state==3) { fEtaBinEdgesv.push_back(edge); }
      else if(state==4) { fDoAbsPhi = (int(edge)==1); state++; }
      else if(state==5) { fPhiBinEdgesv.push_back(edge); }
      else if(state==6) { fNPVBinEdgesv.push_back(edge); }
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
    
  
  unsigned int runNum, lumiSec, evtNum;   // event ID
  unsigned int npv;                       // number of primary vertices
  unsigned int pass;                      // whether probe passes requirements
  float        npu;                       // mean number of expected pileup
  float        scale1fb;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  TLorentzVector *tag=0, *probe=0;        // tag, probe 4-vector
  
  TFile *infile = new TFile(temfname.c_str());   assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
  intree->SetBranchAddress("runNum",   &runNum);
  intree->SetBranchAddress("lumiSec",  &lumiSec);
  intree->SetBranchAddress("evtNum",   &evtNum);
  intree->SetBranchAddress("npv",      &npv);
  intree->SetBranchAddress("npu",      &npu);
  intree->SetBranchAddress("pass",     &pass);
  intree->SetBranchAddress("scale1fb", &scale1fb);
  intree->SetBranchAddress("mass",     &mass);
  intree->SetBranchAddress("qtag",     &qtag);
  intree->SetBranchAddress("qprobe",   &qprobe);
  intree->SetBranchAddress("tag",      &tag);
  intree->SetBranchAddress("probe",    &probe);
  
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
        
        
    if(pass) {
      passPt[ipt]->Fill(mass,puWgt);
      passEta[ieta]->Fill(mass,puWgt);
      passPhi[iphi]->Fill(mass,puWgt);
      passEtaPt[ipt*NBINS_ETA + ieta]->Fill(mass,puWgt);
      passEtaPhi[iphi*NBINS_ETA + ieta]->Fill(mass,puWgt);
      passNPV[inpv]->Fill(mass,puWgt);
    } else {
      failPt[ipt]->Fill(mass,puWgt);
      failEta[ieta]->Fill(mass,puWgt);
      failPhi[iphi]->Fill(mass,puWgt);
      failEtaPt[ipt*NBINS_ETA + ieta]->Fill(mass,puWgt);
      failEtaPhi[iphi*NBINS_ETA + ieta]->Fill(mass,puWgt);
      failNPV[inpv]->Fill(mass,puWgt);
    }    
  }
  infile->Close();
  string outfile_name = fOutputDir + "_binnedTemplates.root";
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
  TLorentzVector *tag=0, *probe=0;        // tag, probe 4-vector
  
  char tname[500];
  
  const unsigned int NBINS_PT  = fPtBinEdgesv.size()-1;
  const unsigned int NBINS_ETA = fEtaBinEdgesv.size()-1;
  const unsigned int NBINS_PHI = fPhiBinEdgesv.size()-1;
  const unsigned int NBINS_NPV = fNPVBinEdgesv.size()-1;
  
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
  
  TFile *infile = new TFile(temfname.c_str());   assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
  intree->SetBranchAddress("runNum",   &runNum);
  intree->SetBranchAddress("lumiSec",  &lumiSec);
  intree->SetBranchAddress("evtNum",   &evtNum);
  intree->SetBranchAddress("npv",      &npv);
  intree->SetBranchAddress("npu",      &npu);
  intree->SetBranchAddress("pass",     &pass);
  intree->SetBranchAddress("scale1fb", &scale1fb);
  intree->SetBranchAddress("mass",     &mass);
  intree->SetBranchAddress("qtag",     &qtag);
  intree->SetBranchAddress("qprobe",   &qprobe);
  intree->SetBranchAddress("tag",      &tag);
  intree->SetBranchAddress("probe",    &probe);  
  
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
        
    if(pass) {
      passPt[ipt]->Fill();
      passEta[ieta]->Fill();
      passPhi[iphi]->Fill();
      passEtaPt[ipt*NBINS_ETA + ieta]->Fill();
      passEtaPhi[iphi*NBINS_ETA + ieta]->Fill();
      passNPV[inpv]->Fill();
    } else {
      failPt[ipt]->Fill();
      failEta[ieta]->Fill();
      failPhi[iphi]->Fill();
      failEtaPt[ipt*NBINS_ETA + ieta]->Fill();
      failEtaPhi[iphi*NBINS_ETA + ieta]->Fill();
      failNPV[inpv]->Fill();
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
    
    if(fSigPass==0 && fSigFail==0) {  // Cut-and-count
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
                   name, cpass, cfail); 
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

      if(fSigPass==0 && fSigFail==0) {  // Cut-and-count
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
                     name, cpass, cfail);
        }
      }

      hEff ->SetBinContent(hEff->GetBin(ix+1, iy+1), eff);
      hErrl->SetBinContent(hErrl->GetBin(ix+1, iy+1), errl);
      hErrh->SetBinContent(hErrh->GetBin(ix+1, iy+1), errh);
    }    
  }  
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
void CEffZFitter::performFit(double &resEff, double &resErrl, double &resErrh,
                             const int ibin,
			     const double xbinLo, const double xbinHi,
			     const double ybinLo, const double ybinHi,
                             TTree *passTree, TTree *failTree,
                             const std::string name, TCanvas *cpass, TCanvas *cfail)
{
  RooRealVar m("m","mass",fFitMassLo,fFitMassHi);
  m.setBins(10000);
  
  char pname[500];
  char binlabelx[1000];
  char binlabely[1000];
  char yield[500];
  char ylabel[500];
  char effstr[1000];
  char nsigstr[1000];
  char nbkgstr[1000];
  char chi2str[1000];
  TFile *histfile = 0;
  if(fSigPass==2 || fSigFail==2) {
    string outfile_name = fOutputDir + "_binnedTemplates.root";
    histfile = new TFile(outfile_name.c_str());
    assert(histfile);
  }
  TFile *datfile = 0;
  if(fSigPass==4 || fSigFail==4) {
    string outfile_name = fOutputDir + "_unbinnedTemplates.root";
    datfile = new TFile(outfile_name.c_str());
    assert(datfile);
  }
  
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
  
  // Define signal and background models
  CSignalModel     *sigModPass = 0;
  CBackgroundModel *bkgModPass = 0;
  CSignalModel     *sigModFail = 0;
  CBackgroundModel *bkgModFail = 0;
  double ptMax=8000;
  double ptMin=0;
  if(name.compare("pt")==0) {
    ptMax=xbinHi;
    ptMin=xbinLo;
  } else if(name.compare("etapt")==0) {
    ptMax=ybinHi;
    ptMin=ybinLo;
  }
 
  if(fSigPass==1) {
    sigModPass = new CBreitWignerConvCrystalBall(m,true);
  
  } else if(fSigPass==2) { 
    char hname[500];
    sprintf(hname,"pass%s_%i",name.c_str(),ibin);
    TH1D *h = (TH1D*)histfile->Get(hname);
    assert(h);
    sigModPass = new CMCTemplateConvGaussian(m,h,true);
  
  } else if(fSigPass==3) {
    sigModPass = new CVoigtianCBShape(m,true);
  
  } else if(fSigPass==4) {
    char tname[500];
    sprintf(tname,"pass%s_%i",name.c_str(),ibin);
    TTree *t = (TTree*)datfile->Get(tname);
    assert(t);
    sigModPass = new CMCDatasetConvGaussian(m,t,true);
  } else if(fSigPass==5) {
    sigModPass = new CBWCBPlusVoigt(m, true, ptMin, ptMax);
  } else if(fSigPass==6) {
    sigModPass = new CBWCBPlusVoigt(m, true, ptMin, ptMax, ibin, name, fSigRefDir);
  }

  if(fBkgPass==1) { 
    bkgModPass = new CExponential(m,true);
  
  } else if(fBkgPass==2) {
    bkgModPass = new CErfcExpo(m,true);
     
  } else if(fBkgPass==3) {
    bkgModPass = new CDoubleExp(m,true);
  
  } else if(fBkgPass==4) {
    bkgModPass = new CLinearExp(m,true);
  
  } else if(fBkgPass==5) {
    bkgModPass = new CQuadraticExp(m,true);
  } else if(fBkgPass==6) {
    bkgModPass = new CErfcExpoFixed(m, true, ibin, name, fBkgRefDir);
  }

  if(fSigFail==1) {
    sigModFail = new CBreitWignerConvCrystalBall(m,false);
  
  } else if(fSigFail==2) {
    char hname[500];
    sprintf(hname,"fail%s_%i",name.c_str(),ibin);
    TH1D *h = (TH1D*)histfile->Get(hname);
    assert(h);
    sigModFail = new CMCTemplateConvGaussian(m,h,false);//,((CMCTemplateConvGaussian*)sigPass)->sigma);

  } else if(fSigFail==3) {
    sigModFail = new CVoigtianCBShape(m,false);
  
  } else if(fSigFail==4) {
    char tname[500];
    sprintf(tname,"fail%s_%i",name.c_str(),ibin);
    TTree *t = (TTree*)datfile->Get(tname);
    assert(t);
    sigModFail = new CMCDatasetConvGaussian(m,t,false);
  } else if(fSigFail==5) {
    sigModFail = new CBWCBPlusVoigt(m, false, ptMin, ptMax);
  } else if(fSigFail==6) {
    sigModFail = new CBWCBPlusVoigt(m, false, ptMin, ptMax, ibin, name, fSigRefDir);
  }
  if(fBkgFail==1) { 
    bkgModFail = new CExponential(m,false);
  
  } else if(fBkgFail==2) {
    bkgModFail = new CErfcExpo(m,false);

  } else if(fBkgFail==3) {
    bkgModFail = new CDoubleExp(m,false);
    
  } else if(fBkgFail==4) {
    bkgModFail = new CLinearExp(m,false);
  
  } else if(fBkgFail==5) {
    bkgModFail = new CQuadraticExp(m,false);
  } else if(fBkgPass==6) {
    bkgModFail = new CErfcExpoFixed(m, false, ibin, name, fBkgRefDir);
  }
  
  // Define free parameters
  double NsigMax     = doBinned ? histPass.Integral()+histFail.Integral() : passTree->GetEntries()+failTree->GetEntries();
  double NbkgFailMax = doBinned ? histFail.Integral() : failTree->GetEntries();
  double NbkgPassMax = doBinned ? histPass.Integral() : passTree->GetEntries();
  // Fine tuning constraints
  //if(name.compare("pt")==0) {
  //  if(xbinLo>=20) NbkgFailMax = NbkgFailMax * 0.3;
  //  else if(xbinLo>=50) NbkgFailMax = NbkgFailMax * 0.05;
  //} else if(name.compare("etapt")==0) {
  //  if(ybinLo>=20) NbkgFailMax = NbkgFailMax * 0.3;
  //  else if(ybinLo>=50) NbkgFailMax = NbkgFailMax * 0.05;
  //}
  RooRealVar Nsig("Nsig","Signal Yield",0.80*NsigMax,0,NsigMax);
  RooRealVar eff("eff","Efficiency",0.8,0,1.0);
  RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",0.01*NbkgPassMax,0,NbkgPassMax);
  RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",0.05*NbkgFailMax,0,NbkgFailMax);  
  if(ptMin<=20) {
    NbkgFail.setVal(0.5*NbkgFailMax);
    //NbkgFail.setMin(0.5*NbkgFailMax);
  }
  if(ptMin>=30) NbkgPass.setMax(NbkgPassMax * 0.02);
  else NbkgPass.setMax(NbkgPassMax * 0.05);
  if(fBkgPass==0) NbkgPass.setVal(0);
  if(fBkgFail==0) NbkgFail.setVal(0);
  RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(eff,Nsig));
  RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));
  RooAddPdf *modelPass=0, *modelFail=0;
  RooExtendPdf *esignalPass=0, *ebackgroundPass=0, *esignalFail=0, *ebackgroundFail=0;
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
                              (fBkgPass>0) ? RooArgList(*(sigModPass->model),*(bkgModPass->model)) :  RooArgList(*(sigModPass->model)),
                              (fBkgPass>0) ? RooArgList(NsigPass,NbkgPass) : RooArgList(NsigPass));
  
    //modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(sigModFail->model),*(bkgModFail->model)),RooArgList(NsigFail,NbkgFail));
    modelFail = new RooAddPdf("modelFail","Model for FAIL sample",
                              (fBkgFail>0) ? RooArgList(*(sigModFail->model),*(bkgModFail->model)) :  RooArgList(*(sigModFail->model)),
                              (fBkgFail>0) ? RooArgList(NsigFail,NbkgFail) : RooArgList(NsigFail));
  }
  
  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");  
  totalPdf.addPdf(*modelFail,"Fail");

  RooFitResult *fitResult=0;
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::Extended(),
                             RooFit::Strategy(2),
                             RooFit::Minos(RooArgSet(eff)),
                             RooFit::NumCPU(4),
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
  
  //
  // Plot passing probes
  //
  sprintf(pname,"pass%s_%i",name.c_str(),ibin);
  sprintf(yield,"%u Events",(int)passTree->GetEntries());
  sprintf(ylabel,"Events / %.1f GeV",(double)BIN_SIZE_PASS);
  sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigPass.getVal(),NsigPass.getPropagatedError(*fitResult));
  sprintf(chi2str,"#chi^{2}/dof = %.3f",mframePass->chiSquare());
  if(fBkgPass>0)
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
  if(fBkgPass>0) {
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
  
  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[1000];    
  sprintf(txtfname,"%s/fitres%s_%i.txt",CPlot::sOutDir.Data(),name.c_str(),ibin);
  printf("\n\nSaved results to %s/fitres%s_%i.txt\n\n",CPlot::sOutDir.Data(),name.c_str(),ibin);
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  fitResult->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResult);
  txtfile.close();
  
  //
  // Clean up
  //
  delete esignalPass;
  delete ebackgroundPass;
  delete esignalFail;
  delete ebackgroundFail;
  delete modelPass;
  delete modelFail;  
  delete dataCombined;
  delete dataPass;
  delete dataFail;
  delete sigModPass;
  delete bkgModPass;
  delete sigModFail;
  delete bkgModFail;        
  delete histfile;
  delete datfile;
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
