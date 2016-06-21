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
  TFile *f_sf_mu        = TFile::Open("~/leptonScaleFactors/root_local/02-02-2016/template_erfcexp/scalefactors_mu.root", "READ");
  TFile *f_data         = TFile::Open("~/leptonScaleFactors/root/SingleMuon_miniAOD_BaselineToDecayModeFinding_tauTnP.root","READ");

  TTree *t_mc_template = (TTree*)f_mc_template->Get("Events");
  TTree *t_data        = (TTree*)f_data->Get("Events");
  
  TH2D *sf_mu_tight = (TH2D*)f_sf_mu->Get("unfactorized_scalefactors_Tight_mu");

  TH1D *h_data_os_pass = new TH1D ("h_data_os_pass", "Opposite sign events in data", 150, 0, 300); 
  TH1D *h_data_os_fail = new TH1D ("h_data_os_fail", "Opposite sign events in data", 150, 0, 300); 
  TH1D *h_data_ss_pass = new TH1D ("h_data_ss_pass", "Same sign events in data", 150, 0, 300); 
  TH1D *h_data_ss_fail = new TH1D ("h_data_ss_fail", "Same sign events in data", 150, 0, 300); 
  TH1D *h_mc_template_pass   = new TH1D ("h_mc_template_pass", "Signal template from truth matched MC", 150,0, 300);
  TH1D *h_mc_template_fail   = new TH1D ("h_mc_template_fail", "Signal template from truth matched MC", 150,0, 300);
  char cuts[512];
  sprintf(cuts, "qtag+qprobe==0 && pass==1 && probe.Pt() > %f && probe.Pt() < %f && TMath::Abs(probe.Eta()) > %f && TMath::Abs(probe.Eta()) < %f", pt_min, pt_max, eta_min, eta_max);
  t_data->Draw("mass>>h_data_os_pass", cuts, "goff");
  sprintf(cuts, "qtag+qprobe==0 && pass==0 && probe.Pt() > %f && probe.Pt() < %f && TMath::Abs(probe.Eta()) > %f && TMath::Abs(probe.Eta()) < %f", pt_min, pt_max, eta_min, eta_max);
  t_data->Draw("mass>>h_data_os_fail", cuts, "goff");
  sprintf(cuts, "qtag+qprobe!=0 && pass==1 && probe.Pt() > %f && probe.Pt() < %f && TMath::Abs(probe.Eta()) > %f && TMath::Abs(probe.Eta()) < %f", pt_min, pt_max, eta_min, eta_max);
  t_data->Draw("mass>>h_data_ss_pass", cuts, "goff");
  sprintf(cuts, "qtag+qprobe!=0 && pass==0 && probe.Pt() > %f && probe.Pt() < %f && TMath::Abs(probe.Eta()) > %f && TMath::Abs(probe.Eta()) < %f", pt_min, pt_max, eta_min, eta_max);
  t_data->Draw("mass>>h_data_ss_fail", cuts, "goff");
 
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
    if(pass==1) h_mc_template_pass->Fill(mass, scale1fb*sf); 
    if(pass==0) h_mc_template_fail->Fill(mass, scale1fb*sf); 
  }

  // do the fit
  TObjArray *shapes_pass = new TObjArray(2);
  shapes_pass->Add(h_data_ss_pass);
  shapes_pass->Add(h_mc_template_pass);
  TFractionFitter *fit_pass = new TFractionFitter(h_data_os_pass, shapes_pass);
  TObjArray *shapes_fail = new TObjArray(2);
  shapes_fail->Add(h_data_ss_fail);
  shapes_fail->Add(h_mc_template_fail);
  TFractionFitter *fit_fail = new TFractionFitter(h_data_os_fail, shapes_fail);

  Int_t status_pass = fit_pass->Fit();
  Int_t status_fail = fit_fail->Fit();
  double coeff_data_ss_pass, d_coeff_data_ss_pass, coeff_data_ss_fail, d_coeff_data_ss_fail;
  double coeff_mc_template_pass, d_coeff_mc_template_pass, coeff_mc_template_fail, d_coeff_mc_template_fail;
  fit_pass->GetResult(1, coeff_data_ss_pass, d_coeff_data_ss_pass);
  fit_fail->GetResult(1, coeff_data_ss_fail, d_coeff_data_ss_fail);
  fit_pass->GetResult(2, coeff_mc_template_pass, d_coeff_mc_template_pass);
  fit_fail->GetResult(2, coeff_mc_template_fail, d_coeff_mc_template_fail);
  if(draw) {
    gStyle->SetOptStat(0);
    char title[512];
    
    h_data_os_pass->SetMarkerStyle(20);
    h_data_os_pass->SetMarkerSize(0.8);
    h_data_os_pass->SetMarkerColor(1);
    h_data_os_pass->SetLineColor(1);
    h_data_ss_pass->SetLineStyle(2);
    h_data_ss_pass->SetLineColor(2);
    h_data_ss_pass->SetLineWidth(2);
    h_mc_template_pass->SetLineStyle(1);
    h_mc_template_pass->SetLineColor(4);
    h_mc_template_pass->SetLineWidth(2);
    TCanvas *c_model_sum_pass = new TCanvas ("c_model_sum_pass", "Passing probes");
    TH1F *model_sum_pass = (TH1F*)fit_pass->GetPlot();
    sprintf(title, "Passing probes, Z #rightarrow #tau_{h}#tau_{#mu} events");
    h_data_os_pass->SetTitle(title);
    h_data_os_pass->Draw("ZP");
    h_data_os_pass->GetXaxis()->SetTitle("mass [GeV]");
    h_data_ss_pass->Draw("L SAME");
    model_sum_pass->SetLineStyle(1);
    model_sum_pass->SetLineColor(4);
    model_sum_pass->SetLineWidth(2);
    model_sum_pass->Draw("L SAME");

    h_data_os_fail->SetMarkerStyle(20);
    h_data_os_fail->SetMarkerSize(0.8);
    h_data_os_fail->SetMarkerColor(1);
    h_data_os_fail->SetLineColor(1);
    h_data_ss_fail->SetLineStyle(2);
    h_data_ss_fail->SetLineColor(2);
    h_data_ss_fail->SetLineWidth(2);
    h_mc_template_fail->SetLineStyle(1);
    h_mc_template_fail->SetLineColor(4);
    h_mc_template_fail->SetLineWidth(2);
    TCanvas *c_model_sum_fail = new TCanvas ("c_model_sum_fail", "Failing probes");
    TH1F *model_sum_fail = (TH1F*)fit_fail->GetPlot();
    sprintf(title, "Failing probes, Z #rightarrow #tau_{h}#tau_{#mu} events");
    h_data_os_fail->SetTitle(title);
    h_data_os_fail->Draw("ZP");
    h_data_os_fail->GetXaxis()->SetTitle("mass [GeV]");
    h_data_ss_fail->Draw("L SAME");
    model_sum_fail->SetLineStyle(1);
    model_sum_fail->SetLineColor(4);
    model_sum_fail->SetLineWidth(2);
    model_sum_fail->Draw("L SAME");
  }
  double efficiency = h_mc_template_pass->Integral() / (h_mc_template_fail->Integral() + h_mc_template_pass->Integral());
  printf("efficiency = %f\n", efficiency);
  return efficiency;
}
