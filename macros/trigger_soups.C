#include <TROOT.h>
#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <vector>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TFile.h>
#include <TCut.h>
#include <TColor.h>
#include <TPaletteAxis.h>
#include <iostream>
#include "leptons.h"
Int_t mit_red  = 1861; 
Int_t mit_gray = 1862; 
Float_t mu_pt_bins[] = {10,20,30,200};
Float_t mu_eta_bins[] = {0., 0.9, 1.2, 2.1, 2.4};
Int_t n_mu_pt_bins=3;
Int_t n_mu_eta_bins=4;
Float_t ele_pt_bins[] = {10,20,30,200};
Float_t ele_eta_bins[] = {0., 1.4442, 1.566, 2.5};
Int_t n_ele_pt_bins=3;
Int_t n_ele_eta_bins=3;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dimuon_dielectron_triggers() {
  TFile *f_dimuon=TFile::Open("~/leptonScaleFactors/root/SingleMuonRun2016B_80x_files_miniAOD_2016-06-01_SoupGivenIsoMu20_soupTnP.root");
  TTree *t_dimuon=(TTree*)f_dimuon->Get("Events");
  TH2D *h_dimuon_eff[9], *h_dimuon_errl[9], *h_dimuon_errh[9];
  char title[512], name[128], cut[512];
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  TPaletteAxis *palette_axis;
  mitPalette();
  for(int i=0; i<n_mu_pt_bins; i++) { for(int j=0; j<n_mu_pt_bins; j++) {
    int nbin = i*n_mu_pt_bins + j;
    printf("Doing pt bin %d ...\n", nbin);
    double tag_pt_min = mu_pt_bins[i],
           tag_pt_max = mu_pt_bins[i+1],
           probe_pt_min = mu_pt_bins[j],
           probe_pt_max = mu_pt_bins[j+1];

    sprintf(name,  "h_dimuon_eff_%d", nbin);
    sprintf(title, "#varepsilon_{soup} (dimuon) tag p_{T} [%d, %d] probe p_{T} [%d, %d]",(int)tag_pt_min, (int)tag_pt_max, (int)probe_pt_min, (int)probe_pt_max);
    h_dimuon_eff[nbin] = new TH2D(name, title, n_mu_eta_bins, mu_eta_bins, n_mu_eta_bins, mu_eta_bins);
    sprintf(name,  "h_dimuon_errl_%d", nbin);
    sprintf(title, "#sigma_{-}^{stat}(#varepsilon_{soup}) (dimuon) tag p_{T} [%d, %d] probe p_{T} [%d, %d]",(int)tag_pt_min, (int)tag_pt_max, (int)probe_pt_min, (int)probe_pt_max);
    h_dimuon_errl[nbin] = new TH2D(name, title, n_mu_eta_bins, mu_eta_bins, n_mu_eta_bins, mu_eta_bins);
    sprintf(name,  "h_dimuon_errh_%d", nbin);
    sprintf(title, "#sigma_{+}^{stat}(#varepsilon_{soup}) (dimuon) tag p_{T} [%d, %d] probe p_{T} [%d, %d]",(int)tag_pt_min, (int)tag_pt_max, (int)probe_pt_min, (int)probe_pt_max);
    h_dimuon_errh[nbin] = new TH2D(name, title, n_mu_eta_bins, mu_eta_bins, n_mu_eta_bins, mu_eta_bins);

    for(int k=0; k<n_mu_eta_bins; k++) { for(int l=0; l<n_mu_eta_bins; l++) {
      double tag_eta_min = mu_eta_bins[k],
             tag_eta_max = mu_eta_bins[k+1],
             probe_eta_min = mu_eta_bins[l],
             probe_eta_max = mu_eta_bins[l+1];
      sprintf(cut, "tag.Pt() < %f && tag.Pt() >= %f && TMath::Abs(tag.Eta()) < %f && TMath::Abs(tag.Eta()) >= %f", tag_pt_max, tag_pt_min, tag_eta_max, tag_eta_min);
      TCut tag_cuts = cut;
      sprintf(cut, "probe.Pt() < %f && probe.Pt() >= %f && TMath::Abs(probe.Eta()) < %f && TMath::Abs(probe.Eta()) >= %f", probe_pt_max, probe_pt_min, probe_eta_max, probe_eta_min);
      TCut probe_cuts = cut;
      TCut passing_requirement = "pass==1";
      TCut passing_events_cut = tag_cuts && probe_cuts && passing_requirement;
      TCut all_events_cut     = tag_cuts && probe_cuts;
      Long64_t passing_events = t_dimuon->GetEntries(passing_events_cut);
      Long64_t all_events     = t_dimuon->GetEntries(all_events_cut);
      double efficiency  = (all_events > 0) ? double(passing_events) / double(all_events) : 0; 
      double error_lo = efficiency - TEfficiency::ClopperPearson(all_events, passing_events, 0.68269, false);
      double error_hi = TEfficiency::ClopperPearson(all_events, passing_events, 0.68269, true) - efficiency;
      int etabin = h_dimuon_eff[nbin]->GetBin(l+1, k+1);
      h_dimuon_eff[nbin]->SetBinContent(etabin, efficiency);
      h_dimuon_eff[nbin]->SetBinError(etabin, TMath::Max(error_lo, error_hi)); //approximate
      h_dimuon_errl[nbin]->SetBinContent(etabin, error_lo);
      h_dimuon_errh[nbin]->SetBinContent(etabin, error_hi);
    }}
    TCanvas *c_dimuon_eff = new TCanvas("c_dimuon_eff", "c_dimuon_eff");
    h_dimuon_eff[nbin]->GetXaxis()->SetTitle("probe |#eta|");
    h_dimuon_eff[nbin]->GetYaxis()->SetTitle("tag |#eta|");
    h_dimuon_eff[nbin]->SetMaximum(1.);
    h_dimuon_eff[nbin]->SetMinimum(0);
    h_dimuon_eff[nbin]->Draw("TEXTE COLZ");
    c_dimuon_eff->Print((string(h_dimuon_eff[nbin]->GetName()) + ".png").c_str());
  }}
}
