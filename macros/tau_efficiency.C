#include <TROOT.h>
#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TFractionFitter.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <vector>
#include <fstream>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLegend.h>
#include <TFile.h>
#include <TCut.h>
#include <leptons.h>
#include <good_run.h>
using namespace std;

double tau_efficiency(double pt_min, double pt_max, double eta_min, double eta_max, bool draw=false) {
  TFile *f_mc_template  = TFile::Open("~/leptonScaleFactors/root/DYJetsToLL_miniAOD_genMatching_BaselineToDecayModeFinding_tauTnP.root", "READ");
  TFile *f_sf_mu        = TFile::Open("~/www/scalefactors/02-02-2016/scalefactors_mu.root", "READ");
  TFile *f_data         = TFile::Open("~/leptonScaleFactors/root/SingleMuon_miniAOD_BaselineToDecayModeFinding_tauTnP.root","READ");

  TTree *t_mc_template = (TTree*)f_mc_template->Get("Events");
  TTree *t_data        = (TTree*)f_data->Get("Events");
  
  TH2D *sf_mu_tight = (TH2D*)f_sf_mu->Get("unfactorized_scalefactors_Tight_mu");

  TH1D *h_data_os = new TH1D ("h_data_os", "Opposite sign events in data", 150, 0, 300); 
  TH1D *h_data_ss = new TH1D ("h_data_ss", "Same sign events in data", 150, 0, 300); 
  TH1D *h_mc_template   = new TH1D ("h_mc_template", "Signal template from truth matched MC", 150,0, 300);
  char cuts[512];
  sprintf(cuts, "qtag+qprobe==0 && probe.Pt() > %f && probe.Pt() < %f && TMath::Abs(probe.Eta()) > %f && TMath::Abs(probe.Eta()) < %f", pt_min, pt_max, eta_min, eta_max);
  t_data->Draw("mass>>h_data_os", cuts, "goff");
  sprintf(cuts, "qtag+qprobe!=0 && probe.Pt() > %f && probe.Pt() < %f && TMath::Abs(probe.Eta()) > %f && TMath::Abs(probe.Eta()) < %f", pt_min, pt_max, eta_min, eta_max);
  t_data->Draw("mass>>h_data_ss", cuts, "goff");
 
  //declare output variables
  unsigned int runNum, // event ID
  lumiSec,
  evtNum,
  npv, // number of primary vertices
  pass; // whether probe passes requirements
  float        npu=1;                     // mean number of expected pileup
  float        scale1fb=1;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  int          truth_tag, truth_probe;              // tag, probe truth
  TLorentzVector *p4_tag=0, *p4_probe=0;        // tag, probe 4-vector 
  t_mc_template->SetBranchAddress("runNum",   &runNum);
  t_mc_template->SetBranchAddress("lumiSec",  &lumiSec);
  t_mc_template->SetBranchAddress("evtNum",   &evtNum);
  t_mc_template->SetBranchAddress("npv",      &npv);
  t_mc_template->SetBranchAddress("pass",     &pass);
  t_mc_template->SetBranchAddress("npu",      &npu);
  t_mc_template->SetBranchAddress("scale1fb", &scale1fb);
  t_mc_template->SetBranchAddress("mass",     &mass);
  t_mc_template->SetBranchAddress("qtag",     &qtag);
  t_mc_template->SetBranchAddress("qprobe",   &qprobe);
  t_mc_template->SetBranchAddress("tag",      &p4_tag);
  t_mc_template->SetBranchAddress("probe",    &p4_probe); 
  t_mc_template->SetBranchAddress("truth_tag",     &truth_tag);
  t_mc_template->SetBranchAddress("truth_probe",   &truth_probe);
  for(Long64_t i = 0; i < t_mc_template->GetEntries(); i++) {
    t_mc_template->GetEntry(i);
    if(!(
      qtag+qprobe==0 &&
      p4_probe->Pt() > pt_min &&
      p4_probe->Pt() < pt_max &&
      TMath::Abs(p4_probe->Eta()) > eta_min &&
      TMath::Abs(p4_probe->Eta()) < eta_max
    )) continue;
    int nbin_sf = sf_mu_tight->FindBin(TMath::Abs(p4_tag->Eta()), p4_tag->Pt());
    double sf = sf_mu_tight->GetBinContent(nbin_sf);
    h_mc_template->Fill(mass, scale1fb*sf); 
  }

  // do the fit
  TObjArray *shapes = new TObjArray(2);
  shapes->Add(h_data_ss);
  shapes->Add(h_mc_template);
  TFractionFitter *fit = new TFractionFitter(h_data_os, shapes);
  Int_t status = fit->Fit();
  if(draw) {
    TCanvas *c_data_os=new TCanvas("c_data_os", "Opposite sign events in data");
    h_data_os->SetMarkerStyle(20);
    h_data_os->SetMarkerSize(0.8);
    h_data_os->SetMarkerColor(1);
    h_data_os->SetLineColor(1);
    h_data_os->Draw("ZP");
    h_data_os->GetXaxis()->SetTitle("mass [GeV]");
    gPad->Update();
    TCanvas *c_data_ss=new TCanvas("c_data_ss", "Same sign events in data");
    h_data_ss->SetLineStyle(2);
    h_data_ss->SetLineColor(2);
    h_data_ss->SetLineWidth(2);
    h_data_ss->Draw("L");
    h_data_ss->GetXaxis()->SetTitle("mass [GeV]");
    gPad->Update();
    TCanvas *c_mc_template=new TCanvas("c_mc_template", "Signal template from truth matched MC");
    h_mc_template->SetLineStyle(1);
    h_mc_template->SetLineColor(4);
    h_mc_template->SetLineWidth(2);
    h_mc_template->Draw("L");
    h_mc_template->GetXaxis()->SetTitle("mass [GeV]");
    gPad->Update();
    TCanvas *c_model_sum = new TCanvas ("c_model_sum", "Model sum");
    TH1F *model_sum = (TH1F*)fit->GetPlot();
    char title[512];
    sprintf(title, "Z #rightarrow #tau_{h}#tau_{#mu} events");
    h_data_os->SetTitle(title);
    h_data_os->Draw("ZP");
    double coeff_data_ss, d_coeff_data_ss;
    fit->GetResult(1, coeff_data_ss, d_coeff_data_ss);
    //h_data_ss->Scale(coeff_data_ss);
    h_data_ss->Draw("L SAME");
    model_sum->SetLineStyle(1);
    model_sum->SetLineColor(4);
    model_sum->SetLineWidth(2);
    model_sum->Draw("L SAME");
  }
  return 0;
}
