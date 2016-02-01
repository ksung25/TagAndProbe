#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TPaletteAxis.h>
#include <iostream>
//Float_t ele_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
//Float_t ele_eta_bins[] = {0, 1.479, 2.4};
//Int_t n_ele_pt_bins=12;
//Int_t n_ele_eta_bins=2;
//Float_t mu_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
//Float_t mu_eta_bins[] = {0, 1.5, 2.1, 2.4};
//Int_t n_mu_pt_bins=12;
//Int_t n_mu_eta_bins=3;
Float_t ele_pt_bins[] = {10,20,30,40,50,70,100,8000};
Float_t ele_eta_bins[] = {0, 1.479, 2.4};
Int_t n_ele_pt_bins=7;
Int_t n_ele_eta_bins=2;
Float_t mu_pt_bins[] = {10,20,30,40,50,70,100,8000};
Float_t mu_eta_bins[] = {0, 1.479, 2.4};
Int_t n_mu_pt_bins=7;
Int_t n_mu_eta_bins=2;

// Function to get MIT colors for the 2D plots
void mitPalette()
{
  static Int_t  colors[100];
  static Bool_t initialized = kFALSE;
  Double_t Red[3]    = { 1, 138./255., 163/255.};
  Double_t Green[3]  = { 1, 139./255., 31/255.};
  Double_t Blue[3]   = { 1, 140./255., 52/255.};
  Double_t Length[3] = { 0.00, 0.35, 1.00 };
  if(!initialized){
    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,100);
    for (int i=0; i<100; i++) colors[i] = FI+i;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(100,colors);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// reads files (names hardcoded) from the root_dir
// takes the two parts of the electron selection and multiplies them
// produces plots and rootfiles for the factorized and unfactorized efficiencies
void factorize_mc_ele(string plots_dir, string root_dir) {
  
  TFile *f_BaselineToLoose_electronTnP          = TFile::Open((plots_dir+"DYJetsToLL_BaselineToLoose_electronTnP/eff.root"       ).c_str() ,"READ");
  TFile *f_BaselineToMedium_electronTnP         = TFile::Open((plots_dir+"DYJetsToLL_BaselineToMedium_electronTnP/eff.root"      ).c_str() ,"READ");
  TFile *f_BaselineToTight_electronTnP          = TFile::Open((plots_dir+"DYJetsToLL_BaselineToTight_electronTnP/eff.root"       ).c_str() ,"READ");
  TFile *f_BaselineToVeto_electronTnP           = TFile::Open((plots_dir+"DYJetsToLL_BaselineToVeto_electronTnP/eff.root"        ).c_str() ,"READ");
  TFile *f_LooseIdToLooseIdIso_electronTnP      = TFile::Open((plots_dir+"DYJetsToLL_LooseIdToLooseIdIso_electronTnP/eff.root"   ).c_str() ,"READ");
  TFile *f_LooseIsoToLooseIdIso_electronTnP     = TFile::Open((plots_dir+"DYJetsToLL_LooseIsoToLooseIdIso_electronTnP/eff.root"  ).c_str() ,"READ");
  TFile *f_MediumIdToMediumIdIso_electronTnP    = TFile::Open((plots_dir+"DYJetsToLL_MediumIdToMediumIdIso_electronTnP/eff.root" ).c_str() ,"READ");
  TFile *f_MediumIsoToMediumIdIso_electronTnP   = TFile::Open((plots_dir+"DYJetsToLL_MediumIsoToMediumIdIso_electronTnP/eff.root").c_str() ,"READ");
  TFile *f_TightIdToTightIdIso_electronTnP      = TFile::Open((plots_dir+"DYJetsToLL_TightIdToTightIdIso_electronTnP/eff.root"   ).c_str() ,"READ");
  TFile *f_TightIsoToTightIdIso_electronTnP     = TFile::Open((plots_dir+"DYJetsToLL_TightIsoToTightIdIso_electronTnP/eff.root"  ).c_str() ,"READ");
  TFile *f_VetoIdToVetoIdIso_electronTnP        = TFile::Open((plots_dir+"DYJetsToLL_VetoIdToVetoIdIso_electronTnP/eff.root"     ).c_str() ,"READ");
  TFile *f_VetoIsoToVetoIdIso_electronTnP       = TFile::Open((plots_dir+"DYJetsToLL_VetoIsoToVetoIdIso_electronTnP/eff.root"    ).c_str() ,"READ");
  
  TH2D * eff_BaselineToLoose_electronTnP        = (TH2D*) f_BaselineToLoose_electronTnP          ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToMedium_electronTnP       = (TH2D*) f_BaselineToMedium_electronTnP         ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToTight_electronTnP        = (TH2D*) f_BaselineToTight_electronTnP          ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToVeto_electronTnP         = (TH2D*) f_BaselineToVeto_electronTnP           ->Get("hEffEtaPt"); 
  TH2D * eff_LooseIdToLooseIdIso_electronTnP    = (TH2D*) f_LooseIdToLooseIdIso_electronTnP      ->Get("hEffEtaPt"); 
  TH2D * eff_LooseIsoToLooseIdIso_electronTnP   = (TH2D*) f_LooseIsoToLooseIdIso_electronTnP     ->Get("hEffEtaPt"); 
  TH2D * eff_MediumIdToMediumIdIso_electronTnP  = (TH2D*) f_MediumIdToMediumIdIso_electronTnP    ->Get("hEffEtaPt"); 
  TH2D * eff_MediumIsoToMediumIdIso_electronTnP = (TH2D*) f_MediumIsoToMediumIdIso_electronTnP   ->Get("hEffEtaPt"); 
  TH2D * eff_TightIdToTightIdIso_electronTnP    = (TH2D*) f_TightIdToTightIdIso_electronTnP      ->Get("hEffEtaPt"); 
  TH2D * eff_TightIsoToTightIdIso_electronTnP   = (TH2D*) f_TightIsoToTightIdIso_electronTnP     ->Get("hEffEtaPt"); 
  TH2D * eff_VetoIdToVetoIdIso_electronTnP      = (TH2D*) f_VetoIdToVetoIdIso_electronTnP        ->Get("hEffEtaPt"); 
  TH2D * eff_VetoIsoToVetoIdIso_electronTnP     = (TH2D*) f_VetoIsoToVetoIdIso_electronTnP       ->Get("hEffEtaPt"); 
  
  TH2D * err_BaselineToLoose_electronTnP        = (TH2D*) f_BaselineToLoose_electronTnP          ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToMedium_electronTnP       = (TH2D*) f_BaselineToMedium_electronTnP         ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToTight_electronTnP        = (TH2D*) f_BaselineToTight_electronTnP          ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToVeto_electronTnP         = (TH2D*) f_BaselineToVeto_electronTnP           ->Get("hErrlEtaPt"); 
  TH2D * err_LooseIdToLooseIdIso_electronTnP    = (TH2D*) f_LooseIdToLooseIdIso_electronTnP      ->Get("hErrlEtaPt"); 
  TH2D * err_LooseIsoToLooseIdIso_electronTnP   = (TH2D*) f_LooseIsoToLooseIdIso_electronTnP     ->Get("hErrlEtaPt"); 
  TH2D * err_MediumIdToMediumIdIso_electronTnP  = (TH2D*) f_MediumIdToMediumIdIso_electronTnP    ->Get("hErrlEtaPt"); 
  TH2D * err_MediumIsoToMediumIdIso_electronTnP = (TH2D*) f_MediumIsoToMediumIdIso_electronTnP   ->Get("hErrlEtaPt"); 
  TH2D * err_TightIdToTightIdIso_electronTnP    = (TH2D*) f_TightIdToTightIdIso_electronTnP      ->Get("hErrlEtaPt"); 
  TH2D * err_TightIsoToTightIdIso_electronTnP   = (TH2D*) f_TightIsoToTightIdIso_electronTnP     ->Get("hErrlEtaPt"); 
  TH2D * err_VetoIdToVetoIdIso_electronTnP      = (TH2D*) f_VetoIdToVetoIdIso_electronTnP        ->Get("hErrlEtaPt"); 
  TH2D * err_VetoIsoToVetoIdIso_electronTnP     = (TH2D*) f_VetoIsoToVetoIdIso_electronTnP       ->Get("hErrlEtaPt"); 

  TH2D * absdiff_Veto_ele   = new TH2D ( "absdiff_Veto_ele",   "absdiff_Veto_ele",   n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * absdiff_Loose_ele  = new TH2D ( "absdiff_Loose_ele",  "absdiff_Loose_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * absdiff_Medium_ele = new TH2D ( "absdiff_Medium_ele", "absdiff_Medium_ele", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * absdiff_Tight_ele  = new TH2D ( "absdiff_Tight_ele",  "absdiff_Tight_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * diffErr_Veto_ele   = new TH2D ( "diffErr_Veto_ele",   "diffErr_Veto_ele",   n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * diffErr_Loose_ele  = new TH2D ( "diffErr_Loose_ele",  "diffErr_Loose_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * diffErr_Medium_ele = new TH2D ( "diffErr_Medium_ele", "diffErr_Medium_ele", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * diffErr_Tight_ele  = new TH2D ( "diffErr_Tight_ele",  "diffErr_Tight_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * factEff_Veto_ele   = new TH2D ( "factEff_Veto_ele",   "factEff_Veto_ele",   n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * factEff_Loose_ele  = new TH2D ( "factEff_Loose_ele",  "factEff_Loose_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * factEff_Medium_ele = new TH2D ( "factEff_Medium_ele", "factEff_Medium_ele", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * factEff_Tight_ele  = new TH2D ( "factEff_Tight_ele",  "factEff_Tight_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);

  for(int i_eta = 1; i_eta <= n_ele_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_ele_pt_bins; i_pt++) {
    int n = absdiff_Veto_ele->GetBin(i_eta, i_pt);
    eff_BaselineToLoose_electronTnP        ->SetBinError(n, err_BaselineToLoose_electronTnP        ->GetBinContent(n) ); 
    eff_BaselineToMedium_electronTnP       ->SetBinError(n, err_BaselineToMedium_electronTnP       ->GetBinContent(n) );
    eff_BaselineToTight_electronTnP        ->SetBinError(n, err_BaselineToTight_electronTnP        ->GetBinContent(n) );
    eff_BaselineToVeto_electronTnP         ->SetBinError(n, err_BaselineToVeto_electronTnP         ->GetBinContent(n) );
    eff_LooseIdToLooseIdIso_electronTnP    ->SetBinError(n, err_LooseIdToLooseIdIso_electronTnP    ->GetBinContent(n) );
    eff_LooseIsoToLooseIdIso_electronTnP   ->SetBinError(n, err_LooseIsoToLooseIdIso_electronTnP   ->GetBinContent(n) );
    eff_MediumIdToMediumIdIso_electronTnP  ->SetBinError(n, err_MediumIdToMediumIdIso_electronTnP  ->GetBinContent(n) );
    eff_MediumIsoToMediumIdIso_electronTnP ->SetBinError(n, err_MediumIsoToMediumIdIso_electronTnP ->GetBinContent(n) );
    eff_TightIdToTightIdIso_electronTnP    ->SetBinError(n, err_TightIdToTightIdIso_electronTnP    ->GetBinContent(n) );
    eff_TightIsoToTightIdIso_electronTnP   ->SetBinError(n, err_TightIsoToTightIdIso_electronTnP   ->GetBinContent(n) );
    eff_VetoIdToVetoIdIso_electronTnP      ->SetBinError(n, err_VetoIdToVetoIdIso_electronTnP      ->GetBinContent(n) );
    eff_VetoIsoToVetoIdIso_electronTnP     ->SetBinError(n, err_VetoIsoToVetoIdIso_electronTnP     ->GetBinContent(n) );
  }}
  
  (*factEff_Veto_ele  ) = (*eff_VetoIdToVetoIdIso_electronTnP     )  * (*eff_VetoIsoToVetoIdIso_electronTnP      );
  (*factEff_Loose_ele ) = (*eff_LooseIdToLooseIdIso_electronTnP   )  * (*eff_LooseIsoToLooseIdIso_electronTnP    );
  (*factEff_Medium_ele) = (*eff_MediumIdToMediumIdIso_electronTnP )  * (*eff_MediumIsoToMediumIdIso_electronTnP  );
  (*factEff_Tight_ele ) = (*eff_TightIdToTightIdIso_electronTnP   )  * (*eff_TightIsoToTightIdIso_electronTnP    );
 
  (*absdiff_Veto_ele  ) =  (*factEff_Veto_ele   ) - (*eff_BaselineToVeto_electronTnP    );
  (*absdiff_Loose_ele ) =  (*factEff_Loose_ele  ) - (*eff_BaselineToLoose_electronTnP   );
  (*absdiff_Medium_ele) =  (*factEff_Medium_ele ) - (*eff_BaselineToMedium_electronTnP  );
  (*absdiff_Tight_ele ) =  (*factEff_Tight_ele  ) - (*eff_BaselineToTight_electronTnP   );

  // take absolute value of difference bins, and compute the statistical error
  // electrons first
  for(int i_eta = 1; i_eta <= n_ele_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_ele_pt_bins; i_pt++) {
    
    int n = absdiff_Veto_ele->GetBin(i_eta, i_pt);
    
    absdiff_Veto_ele   ->SetBinContent(n, fabs(absdiff_Veto_ele   ->GetBinContent(n) ));
    absdiff_Loose_ele  ->SetBinContent(n, fabs(absdiff_Loose_ele  ->GetBinContent(n) ));
    absdiff_Medium_ele ->SetBinContent(n, fabs(absdiff_Medium_ele ->GetBinContent(n) ));
    absdiff_Tight_ele  ->SetBinContent(n, fabs(absdiff_Tight_ele  ->GetBinContent(n) ));
    
    diffErr_Veto_ele->SetBinContent(n, sqrt(
      pow(err_BaselineToVeto_electronTnP->GetBinContent(n), 2) +
      pow(eff_VetoIsoToVetoIdIso_electronTnP->GetBinContent(n) * err_VetoIdToVetoIdIso_electronTnP->GetBinContent(n),2) +
      pow(eff_VetoIdToVetoIdIso_electronTnP->GetBinContent(n) * err_VetoIsoToVetoIdIso_electronTnP->GetBinContent(n),2)
    ));
    diffErr_Loose_ele->SetBinContent(n, sqrt(
      pow(err_BaselineToLoose_electronTnP->GetBinContent(n), 2) +
      pow(eff_LooseIsoToLooseIdIso_electronTnP->GetBinContent(n) * err_LooseIdToLooseIdIso_electronTnP->GetBinContent(n),2) +
      pow(eff_LooseIdToLooseIdIso_electronTnP->GetBinContent(n) * err_LooseIsoToLooseIdIso_electronTnP->GetBinContent(n),2)
    ));
    diffErr_Medium_ele->SetBinContent(n, sqrt(
      pow(err_BaselineToMedium_electronTnP->GetBinContent(n), 2) +
      pow(eff_MediumIsoToMediumIdIso_electronTnP->GetBinContent(n) * err_MediumIdToMediumIdIso_electronTnP->GetBinContent(n),2) +
      pow(eff_MediumIdToMediumIdIso_electronTnP->GetBinContent(n) * err_MediumIsoToMediumIdIso_electronTnP->GetBinContent(n),2)
    ));
    diffErr_Tight_ele->SetBinContent(n, sqrt(
      pow(err_BaselineToTight_electronTnP->GetBinContent(n), 2) +
      pow(eff_TightIsoToTightIdIso_electronTnP->GetBinContent(n) * err_TightIdToTightIdIso_electronTnP->GetBinContent(n),2) +
      pow(eff_TightIdToTightIdIso_electronTnP->GetBinContent(n) * err_TightIsoToTightIdIso_electronTnP->GetBinContent(n),2)
    ));
  }}

  // Start drawing stuff 
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;

  // draw |difference| plots
  TCanvas *c_absdiff_Veto_ele = new TCanvas("c_absdiff_Veto_ele","Difference from factorized efficiency, Veto electrons (MC)",800,800);
  absdiff_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Veto_ele->SetTitle("Difference from factorized efficiency, Veto electrons (MC)");
  absdiff_Veto_ele->GetXaxis()->SetTitle("| #eta |");
  absdiff_Veto_ele->GetXaxis()->SetTitleOffset(0.9);
  absdiff_Veto_ele->GetXaxis()->SetTitleSize(0.04);
  absdiff_Veto_ele->GetXaxis()->SetLabelSize(0.02);
  absdiff_Veto_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absdiff_Veto_ele->GetYaxis()->SetTitleOffset(0.9);
  absdiff_Veto_ele->GetYaxis()->SetTitleSize(0.04);
  absdiff_Veto_ele->GetYaxis()->SetLabelSize(0.02);
  absdiff_Veto_ele->GetYaxis()->SetRangeUser(10,100);
  absdiff_Veto_ele->SetMinimum(0);
  absdiff_Veto_ele->SetMaximum(0.2);
  absdiff_Veto_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absdiff_Veto_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absdiff_Veto_ele->Update();
  c_absdiff_Veto_ele->Print((plots_dir+"DYJetsToLL_absdiff_Veto_ele.png").c_str());
  TCanvas *c_absdiff_Loose_ele = new TCanvas("c_absdiff_Loose_ele","Difference from factorized efficiency, Loose electrons (MC)",800,800);
  absdiff_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Loose_ele->SetTitle("Difference from factorized efficiency, Loose electrons (MC)");
  absdiff_Loose_ele->GetXaxis()->SetTitle("| #eta |");
  absdiff_Loose_ele->GetXaxis()->SetTitleOffset(0.9);
  absdiff_Loose_ele->GetXaxis()->SetTitleSize(0.04);
  absdiff_Loose_ele->GetXaxis()->SetLabelSize(0.02);
  absdiff_Loose_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absdiff_Loose_ele->GetYaxis()->SetTitleOffset(0.9);
  absdiff_Loose_ele->GetYaxis()->SetTitleSize(0.04);
  absdiff_Loose_ele->GetYaxis()->SetLabelSize(0.02);
  absdiff_Loose_ele->GetYaxis()->SetRangeUser(10,100);
  absdiff_Loose_ele->SetMinimum(0);
  absdiff_Loose_ele->SetMaximum(0.2);
  absdiff_Loose_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absdiff_Loose_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absdiff_Loose_ele->Update();
  c_absdiff_Loose_ele->Print((plots_dir+"DYJetsToLL_absdiff_Loose_ele.png").c_str());
  TCanvas *c_absdiff_Medium_ele = new TCanvas("c_absdiff_Medium_ele","Difference from factorized efficiency, Medium electrons (MC)",800,800);
  absdiff_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Medium_ele->SetTitle("Difference from factorized efficiency, Medium electrons ID (MC)");
  absdiff_Medium_ele->GetXaxis()->SetTitle("| #eta |");
  absdiff_Medium_ele->GetXaxis()->SetTitleOffset(0.9);
  absdiff_Medium_ele->GetXaxis()->SetTitleSize(0.04);
  absdiff_Medium_ele->GetXaxis()->SetLabelSize(0.02);
  absdiff_Medium_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absdiff_Medium_ele->GetYaxis()->SetTitleOffset(0.9);
  absdiff_Medium_ele->GetYaxis()->SetTitleSize(0.04);
  absdiff_Medium_ele->GetYaxis()->SetLabelSize(0.02);
  absdiff_Medium_ele->GetYaxis()->SetRangeUser(10,100);
  absdiff_Medium_ele->SetMinimum(0);
  absdiff_Medium_ele->SetMaximum(0.2);
  absdiff_Medium_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absdiff_Medium_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absdiff_Medium_ele->Update();
  c_absdiff_Medium_ele->Print((plots_dir+"DYJetsToLL_absdiff_Medium_ele.png").c_str());
  TCanvas *c_absdiff_Tight_ele = new TCanvas("c_absdiff_Tight_ele","Difference from factorized efficiency, Tight electrons (MC)",800,800);
  absdiff_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Tight_ele->SetTitle("Difference from factorized efficiency, Tight electrons (MC)");
  absdiff_Tight_ele->GetXaxis()->SetTitle("| #eta |");
  absdiff_Tight_ele->GetXaxis()->SetTitleOffset(0.9);
  absdiff_Tight_ele->GetXaxis()->SetTitleSize(0.04);
  absdiff_Tight_ele->GetXaxis()->SetLabelSize(0.02);
  absdiff_Tight_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absdiff_Tight_ele->GetYaxis()->SetTitleOffset(0.9);
  absdiff_Tight_ele->GetYaxis()->SetTitleSize(0.04);
  absdiff_Tight_ele->GetYaxis()->SetLabelSize(0.02);
  absdiff_Tight_ele->GetYaxis()->SetRangeUser(10,100);
  absdiff_Tight_ele->SetMinimum(0);
  absdiff_Tight_ele->SetMaximum(0.2);
  absdiff_Tight_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absdiff_Tight_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absdiff_Tight_ele->Update();
  c_absdiff_Tight_ele->Print((plots_dir+"DYJetsToLL_absdiff_Tight_ele.png").c_str());

  // Draw unfactorized efficiency plots
  TCanvas *c_eff_Veto_ele = new TCanvas("c_eff_Veto_ele","Unfactorized Veto electron efficiency (MC)",800,800);
  eff_BaselineToVeto_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToVeto_electronTnP->SetTitle("Unfactorized Veto electron efficiency (MC)");
  eff_BaselineToVeto_electronTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToVeto_electronTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToVeto_electronTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToVeto_electronTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToVeto_electronTnP->SetMinimum(0.0);
  eff_BaselineToVeto_electronTnP->SetMaximum(1);
  eff_BaselineToVeto_electronTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToVeto_electronTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Veto_ele->Update();
  c_eff_Veto_ele->Print((plots_dir+"DYJetsToLL_eff_Veto_ele.png").c_str());
  TCanvas *c_eff_Loose_ele = new TCanvas("c_eff_Loose_ele","Unfactorized Loose electron efficiency (MC)",800,800);
  eff_BaselineToLoose_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToLoose_electronTnP->SetTitle("Unfactorized Loose electron efficiency (MC)");
  eff_BaselineToLoose_electronTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToLoose_electronTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToLoose_electronTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToLoose_electronTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToLoose_electronTnP->SetMinimum(0.0);
  eff_BaselineToLoose_electronTnP->SetMaximum(1);
  eff_BaselineToLoose_electronTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToLoose_electronTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Loose_ele->Update();
  c_eff_Loose_ele->Print((plots_dir+"DYJetsToLL_eff_Loose_ele.png").c_str());
  TCanvas *c_eff_Medium_ele = new TCanvas("c_eff_Medium_ele","Unfactorized Medium electron efficiency (MC)",800,800);
  eff_BaselineToMedium_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToMedium_electronTnP->SetTitle("Unfactorized Medium electron efficiency (MC)");
  eff_BaselineToMedium_electronTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToMedium_electronTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToMedium_electronTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToMedium_electronTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToMedium_electronTnP->SetMinimum(0.0);
  eff_BaselineToMedium_electronTnP->SetMaximum(1);
  eff_BaselineToMedium_electronTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToMedium_electronTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Medium_ele->Update();
  c_eff_Medium_ele->Print((plots_dir+"DYJetsToLL_eff_Medium_ele.png").c_str());
  TCanvas *c_eff_Tight_ele = new TCanvas("c_eff_Tight_ele","Unfactorized Tight electron efficiency (MC)",800,800);
  eff_BaselineToTight_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToTight_electronTnP->SetTitle("Unfactorized Tight electron efficiency (MC)");
  eff_BaselineToTight_electronTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToTight_electronTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToTight_electronTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToTight_electronTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToTight_electronTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToTight_electronTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToTight_electronTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToTight_electronTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToTight_electronTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToTight_electronTnP->SetMinimum(0.0);
  eff_BaselineToTight_electronTnP->SetMaximum(1);
  eff_BaselineToTight_electronTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToTight_electronTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Tight_ele->Update();
  c_eff_Tight_ele->Print((plots_dir+"DYJetsToLL_eff_Tight_ele.png").c_str());

  // Draw factorized efficiency plots
  TCanvas *c_factEff_Veto_ele = new TCanvas("c_factEff_Veto_ele","Factorized Veto electron efficiency (MC)",800,800);
  factEff_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Veto_ele->SetTitle("Factorized Veto electron efficiency (MC)");
  factEff_Veto_ele->GetXaxis()->SetTitle("| #eta |");
  factEff_Veto_ele->GetXaxis()->SetTitleOffset(0.9);
  factEff_Veto_ele->GetXaxis()->SetTitleSize(0.04);
  factEff_Veto_ele->GetXaxis()->SetLabelSize(0.02);
  factEff_Veto_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factEff_Veto_ele->GetYaxis()->SetTitleOffset(0.9);
  factEff_Veto_ele->GetYaxis()->SetTitleSize(0.04);
  factEff_Veto_ele->GetYaxis()->SetLabelSize(0.02);
  factEff_Veto_ele->GetYaxis()->SetRangeUser(10,100);
  factEff_Veto_ele->SetMinimum(0.0);
  factEff_Veto_ele->SetMaximum(1);
  factEff_Veto_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factEff_Veto_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factEff_Veto_ele->Update();
  c_factEff_Veto_ele->Print((plots_dir+"DYJetsToLL_factEff_Veto_ele.png").c_str());
  TCanvas *c_factEff_Loose_ele = new TCanvas("c_factEff_Loose_ele","Factorized Loose electron efficiency (MC)",800,800);
  factEff_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Loose_ele->SetTitle("Factorized Loose electron efficiency (MC)");
  factEff_Loose_ele->GetXaxis()->SetTitle("| #eta |");
  factEff_Loose_ele->GetXaxis()->SetTitleOffset(0.9);
  factEff_Loose_ele->GetXaxis()->SetTitleSize(0.04);
  factEff_Loose_ele->GetXaxis()->SetLabelSize(0.02);
  factEff_Loose_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factEff_Loose_ele->GetYaxis()->SetTitleOffset(0.9);
  factEff_Loose_ele->GetYaxis()->SetTitleSize(0.04);
  factEff_Loose_ele->GetYaxis()->SetLabelSize(0.02);
  factEff_Loose_ele->GetYaxis()->SetRangeUser(10,100);
  factEff_Loose_ele->SetMinimum(0.0);
  factEff_Loose_ele->SetMaximum(1);
  factEff_Loose_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factEff_Loose_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factEff_Loose_ele->Update();
  c_factEff_Loose_ele->Print((plots_dir+"DYJetsToLL_factEff_Loose_ele.png").c_str());
  TCanvas *c_factEff_Medium_ele = new TCanvas("c_factEff_Medium_ele","Factorized Medium electron efficiency (MC)",800,800);
  factEff_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Medium_ele->SetTitle("Factorized Medium electron efficiency (MC)");
  factEff_Medium_ele->GetXaxis()->SetTitle("| #eta |");
  factEff_Medium_ele->GetXaxis()->SetTitleOffset(0.9);
  factEff_Medium_ele->GetXaxis()->SetTitleSize(0.04);
  factEff_Medium_ele->GetXaxis()->SetLabelSize(0.02);
  factEff_Medium_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factEff_Medium_ele->GetYaxis()->SetTitleOffset(0.9);
  factEff_Medium_ele->GetYaxis()->SetTitleSize(0.04);
  factEff_Medium_ele->GetYaxis()->SetLabelSize(0.02);
  factEff_Medium_ele->GetYaxis()->SetRangeUser(10,100);
  factEff_Medium_ele->SetMinimum(0.0);
  factEff_Medium_ele->SetMaximum(1);
  factEff_Medium_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factEff_Medium_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factEff_Medium_ele->Update();
  c_factEff_Medium_ele->Print((plots_dir+"DYJetsToLL_factEff_Medium_ele.png").c_str());
  TCanvas *c_factEff_Tight_ele = new TCanvas("c_factEff_Tight_ele","Factorized Tight electron efficiency (MC)",800,800);
  factEff_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Tight_ele->SetTitle("Factorized Tight electron efficiency (MC)");
  factEff_Tight_ele->GetXaxis()->SetTitle("| #eta |");
  factEff_Tight_ele->GetXaxis()->SetTitleOffset(0.9);
  factEff_Tight_ele->GetXaxis()->SetTitleSize(0.04);
  factEff_Tight_ele->GetXaxis()->SetLabelSize(0.02);
  factEff_Tight_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factEff_Tight_ele->GetYaxis()->SetTitleOffset(0.9);
  factEff_Tight_ele->GetYaxis()->SetTitleSize(0.04);
  factEff_Tight_ele->GetYaxis()->SetLabelSize(0.02);
  factEff_Tight_ele->GetYaxis()->SetRangeUser(10,100);
  factEff_Tight_ele->SetMinimum(0.0);
  factEff_Tight_ele->SetMaximum(1);
  factEff_Tight_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factEff_Tight_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factEff_Tight_ele->Update();
  c_factEff_Tight_ele->Print((plots_dir+"DYJetsToLL_factEff_Tight_ele.png").c_str());
  
  TFile *efficiencies_file = TFile::Open((root_dir+"DYJetsToLL_efficiencies_electronTnP.root").c_str(),"RECREATE");
  eff_BaselineToVeto_electronTnP   ->Write("eff_Veto_ele");  
  eff_BaselineToLoose_electronTnP  ->Write("eff_Loose_ele"); 
  eff_BaselineToMedium_electronTnP ->Write("eff_Medium_ele"); 
  eff_BaselineToTight_electronTnP  ->Write("eff_Tight_ele");  
  factEff_Veto_ele   ->Write("factEff_Veto_ele");
  factEff_Loose_ele  ->Write("factEff_Loose_ele");
  factEff_Medium_ele ->Write("factEff_Medium_ele");
  factEff_Tight_ele  ->Write("factEff_Tight_ele"); 
  efficiencies_file->Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void factorize_data_ele(string plots_dir, string root_dir) {
  
  TFile *f_BaselineToLoose_electronTnP          = TFile::Open((plots_dir+"SingleElectron_BaselineToLoose_electronTnP/eff.root"       ).c_str() ,"READ");
  TFile *f_BaselineToMedium_electronTnP         = TFile::Open((plots_dir+"SingleElectron_BaselineToMedium_electronTnP/eff.root"      ).c_str() ,"READ");
  TFile *f_BaselineToTight_electronTnP          = TFile::Open((plots_dir+"SingleElectron_BaselineToTight_electronTnP/eff.root"       ).c_str() ,"READ");
  TFile *f_BaselineToVeto_electronTnP           = TFile::Open((plots_dir+"SingleElectron_BaselineToVeto_electronTnP/eff.root"        ).c_str() ,"READ");
  TFile *f_LooseIdToLooseIdIso_electronTnP      = TFile::Open((plots_dir+"SingleElectron_LooseIdToLooseIdIso_electronTnP/eff.root"   ).c_str() ,"READ");
  TFile *f_LooseIsoToLooseIdIso_electronTnP     = TFile::Open((plots_dir+"SingleElectron_LooseIsoToLooseIdIso_electronTnP/eff.root"  ).c_str() ,"READ");
  TFile *f_MediumIdToMediumIdIso_electronTnP    = TFile::Open((plots_dir+"SingleElectron_MediumIdToMediumIdIso_electronTnP/eff.root" ).c_str() ,"READ");
  TFile *f_MediumIsoToMediumIdIso_electronTnP   = TFile::Open((plots_dir+"SingleElectron_MediumIsoToMediumIdIso_electronTnP/eff.root").c_str() ,"READ");
  TFile *f_TightIdToTightIdIso_electronTnP      = TFile::Open((plots_dir+"SingleElectron_TightIdToTightIdIso_electronTnP/eff.root"   ).c_str() ,"READ");
  TFile *f_TightIsoToTightIdIso_electronTnP     = TFile::Open((plots_dir+"SingleElectron_TightIsoToTightIdIso_electronTnP/eff.root"  ).c_str() ,"READ");
  TFile *f_VetoIdToVetoIdIso_electronTnP        = TFile::Open((plots_dir+"SingleElectron_VetoIdToVetoIdIso_electronTnP/eff.root"     ).c_str() ,"READ");
  TFile *f_VetoIsoToVetoIdIso_electronTnP       = TFile::Open((plots_dir+"SingleElectron_VetoIsoToVetoIdIso_electronTnP/eff.root"    ).c_str() ,"READ");
  
  TH2D * eff_BaselineToLoose_electronTnP        = (TH2D*) f_BaselineToLoose_electronTnP          ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToMedium_electronTnP       = (TH2D*) f_BaselineToMedium_electronTnP         ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToTight_electronTnP        = (TH2D*) f_BaselineToTight_electronTnP          ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToVeto_electronTnP         = (TH2D*) f_BaselineToVeto_electronTnP           ->Get("hEffEtaPt"); 
  TH2D * eff_LooseIdToLooseIdIso_electronTnP    = (TH2D*) f_LooseIdToLooseIdIso_electronTnP      ->Get("hEffEtaPt"); 
  TH2D * eff_LooseIsoToLooseIdIso_electronTnP   = (TH2D*) f_LooseIsoToLooseIdIso_electronTnP     ->Get("hEffEtaPt"); 
  TH2D * eff_MediumIdToMediumIdIso_electronTnP  = (TH2D*) f_MediumIdToMediumIdIso_electronTnP    ->Get("hEffEtaPt"); 
  TH2D * eff_MediumIsoToMediumIdIso_electronTnP = (TH2D*) f_MediumIsoToMediumIdIso_electronTnP   ->Get("hEffEtaPt"); 
  TH2D * eff_TightIdToTightIdIso_electronTnP    = (TH2D*) f_TightIdToTightIdIso_electronTnP      ->Get("hEffEtaPt"); 
  TH2D * eff_TightIsoToTightIdIso_electronTnP   = (TH2D*) f_TightIsoToTightIdIso_electronTnP     ->Get("hEffEtaPt"); 
  TH2D * eff_VetoIdToVetoIdIso_electronTnP      = (TH2D*) f_VetoIdToVetoIdIso_electronTnP        ->Get("hEffEtaPt"); 
  TH2D * eff_VetoIsoToVetoIdIso_electronTnP     = (TH2D*) f_VetoIsoToVetoIdIso_electronTnP       ->Get("hEffEtaPt"); 
  
  TH2D * err_BaselineToLoose_electronTnP        = (TH2D*) f_BaselineToLoose_electronTnP          ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToMedium_electronTnP       = (TH2D*) f_BaselineToMedium_electronTnP         ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToTight_electronTnP        = (TH2D*) f_BaselineToTight_electronTnP          ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToVeto_electronTnP         = (TH2D*) f_BaselineToVeto_electronTnP           ->Get("hErrlEtaPt"); 
  TH2D * err_LooseIdToLooseIdIso_electronTnP    = (TH2D*) f_LooseIdToLooseIdIso_electronTnP      ->Get("hErrlEtaPt"); 
  TH2D * err_LooseIsoToLooseIdIso_electronTnP   = (TH2D*) f_LooseIsoToLooseIdIso_electronTnP     ->Get("hErrlEtaPt"); 
  TH2D * err_MediumIdToMediumIdIso_electronTnP  = (TH2D*) f_MediumIdToMediumIdIso_electronTnP    ->Get("hErrlEtaPt"); 
  TH2D * err_MediumIsoToMediumIdIso_electronTnP = (TH2D*) f_MediumIsoToMediumIdIso_electronTnP   ->Get("hErrlEtaPt"); 
  TH2D * err_TightIdToTightIdIso_electronTnP    = (TH2D*) f_TightIdToTightIdIso_electronTnP      ->Get("hErrlEtaPt"); 
  TH2D * err_TightIsoToTightIdIso_electronTnP   = (TH2D*) f_TightIsoToTightIdIso_electronTnP     ->Get("hErrlEtaPt"); 
  TH2D * err_VetoIdToVetoIdIso_electronTnP      = (TH2D*) f_VetoIdToVetoIdIso_electronTnP        ->Get("hErrlEtaPt"); 
  TH2D * err_VetoIsoToVetoIdIso_electronTnP     = (TH2D*) f_VetoIsoToVetoIdIso_electronTnP       ->Get("hErrlEtaPt"); 

  TH2D * absdiff_Veto_ele   = new TH2D ( "absdiff_Veto_ele",   "absdiff_Veto_ele",   n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * absdiff_Loose_ele  = new TH2D ( "absdiff_Loose_ele",  "absdiff_Loose_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * absdiff_Medium_ele = new TH2D ( "absdiff_Medium_ele", "absdiff_Medium_ele", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * absdiff_Tight_ele  = new TH2D ( "absdiff_Tight_ele",  "absdiff_Tight_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * diffErr_Veto_ele   = new TH2D ( "diffErr_Veto_ele",   "diffErr_Veto_ele",   n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * diffErr_Loose_ele  = new TH2D ( "diffErr_Loose_ele",  "diffErr_Loose_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * diffErr_Medium_ele = new TH2D ( "diffErr_Medium_ele", "diffErr_Medium_ele", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * diffErr_Tight_ele  = new TH2D ( "diffErr_Tight_ele",  "diffErr_Tight_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * factEff_Veto_ele   = new TH2D ( "factEff_Veto_ele",   "factEff_Veto_ele",   n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * factEff_Loose_ele  = new TH2D ( "factEff_Loose_ele",  "factEff_Loose_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * factEff_Medium_ele = new TH2D ( "factEff_Medium_ele", "factEff_Medium_ele", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D * factEff_Tight_ele  = new TH2D ( "factEff_Tight_ele",  "factEff_Tight_ele",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);

  for(int i_eta = 1; i_eta <= n_ele_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_ele_pt_bins; i_pt++) {
    int n = absdiff_Veto_ele->GetBin(i_eta, i_pt);
    eff_BaselineToLoose_electronTnP        ->SetBinError(n, err_BaselineToLoose_electronTnP        ->GetBinContent(n) ); 
    eff_BaselineToMedium_electronTnP       ->SetBinError(n, err_BaselineToMedium_electronTnP       ->GetBinContent(n) );
    eff_BaselineToTight_electronTnP        ->SetBinError(n, err_BaselineToTight_electronTnP        ->GetBinContent(n) );
    eff_BaselineToVeto_electronTnP         ->SetBinError(n, err_BaselineToVeto_electronTnP         ->GetBinContent(n) );
    eff_LooseIdToLooseIdIso_electronTnP    ->SetBinError(n, err_LooseIdToLooseIdIso_electronTnP    ->GetBinContent(n) );
    eff_LooseIsoToLooseIdIso_electronTnP   ->SetBinError(n, err_LooseIsoToLooseIdIso_electronTnP   ->GetBinContent(n) );
    eff_MediumIdToMediumIdIso_electronTnP  ->SetBinError(n, err_MediumIdToMediumIdIso_electronTnP  ->GetBinContent(n) );
    eff_MediumIsoToMediumIdIso_electronTnP ->SetBinError(n, err_MediumIsoToMediumIdIso_electronTnP ->GetBinContent(n) );
    eff_TightIdToTightIdIso_electronTnP    ->SetBinError(n, err_TightIdToTightIdIso_electronTnP    ->GetBinContent(n) );
    eff_TightIsoToTightIdIso_electronTnP   ->SetBinError(n, err_TightIsoToTightIdIso_electronTnP   ->GetBinContent(n) );
    eff_VetoIdToVetoIdIso_electronTnP      ->SetBinError(n, err_VetoIdToVetoIdIso_electronTnP      ->GetBinContent(n) );
    eff_VetoIsoToVetoIdIso_electronTnP     ->SetBinError(n, err_VetoIsoToVetoIdIso_electronTnP     ->GetBinContent(n) );
  }}
  
  (*factEff_Veto_ele  ) = (*eff_VetoIdToVetoIdIso_electronTnP     )  * (*eff_VetoIsoToVetoIdIso_electronTnP      );
  (*factEff_Loose_ele ) = (*eff_LooseIdToLooseIdIso_electronTnP   )  * (*eff_LooseIsoToLooseIdIso_electronTnP    );
  (*factEff_Medium_ele) = (*eff_MediumIdToMediumIdIso_electronTnP )  * (*eff_MediumIsoToMediumIdIso_electronTnP  );
  (*factEff_Tight_ele ) = (*eff_TightIdToTightIdIso_electronTnP   )  * (*eff_TightIsoToTightIdIso_electronTnP    );
 
  // sanity check of error propagation
  //printf("eff_VetoIdToVetoIdIso_electronTnP(1,1) = %f +/- %f\n", eff_VetoIdToVetoIdIso_electronTnP->GetBinContent(eff_VetoIdToVetoIdIso_electronTnP->GetBin(1,1)), eff_VetoIdToVetoIdIso_electronTnP->GetBinError(eff_VetoIdToVetoIdIso_electronTnP->GetBin(1,1)));
  //printf("eff_VetoIsoToVetoIdIso_electronTnP(1,1) = %f +/- %f\n", eff_VetoIsoToVetoIdIso_electronTnP->GetBinContent(eff_VetoIsoToVetoIdIso_electronTnP->GetBin(1,1)), eff_VetoIsoToVetoIdIso_electronTnP->GetBinError(eff_VetoIsoToVetoIdIso_electronTnP->GetBin(1,1)));
  //printf("factEff_Veto_ele(1,1) = %f +/- %f\n", factEff_Veto_ele->GetBinContent(factEff_Veto_ele->GetBin(1,1)), factEff_Veto_ele->GetBinError(factEff_Veto_ele->GetBin(1,1)));
  
  (*absdiff_Veto_ele  ) =  (*factEff_Veto_ele   ) - (*eff_BaselineToVeto_electronTnP    );
  (*absdiff_Loose_ele ) =  (*factEff_Loose_ele  ) - (*eff_BaselineToLoose_electronTnP   );
  (*absdiff_Medium_ele) =  (*factEff_Medium_ele ) - (*eff_BaselineToMedium_electronTnP  );
  (*absdiff_Tight_ele ) =  (*factEff_Tight_ele  ) - (*eff_BaselineToTight_electronTnP   );

  // take absolute value of difference bins, and compute the statistical error
  // electrons first
  for(int i_eta = 1; i_eta <= n_ele_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_ele_pt_bins; i_pt++) {
    
    int n = absdiff_Veto_ele->GetBin(i_eta, i_pt);
    
    absdiff_Veto_ele   ->SetBinContent(n, fabs(absdiff_Veto_ele   ->GetBinContent(n) ));
    absdiff_Loose_ele  ->SetBinContent(n, fabs(absdiff_Loose_ele  ->GetBinContent(n) ));
    absdiff_Medium_ele ->SetBinContent(n, fabs(absdiff_Medium_ele ->GetBinContent(n) ));
    absdiff_Tight_ele  ->SetBinContent(n, fabs(absdiff_Tight_ele  ->GetBinContent(n) ));
    
    diffErr_Veto_ele->SetBinContent(n, sqrt(
      pow(err_BaselineToVeto_electronTnP->GetBinContent(n), 2) +
      pow(eff_VetoIsoToVetoIdIso_electronTnP->GetBinContent(n) * err_VetoIdToVetoIdIso_electronTnP->GetBinContent(n),2) +
      pow(eff_VetoIdToVetoIdIso_electronTnP->GetBinContent(n) * err_VetoIsoToVetoIdIso_electronTnP->GetBinContent(n),2)
    ));
    diffErr_Loose_ele->SetBinContent(n, sqrt(
      pow(err_BaselineToLoose_electronTnP->GetBinContent(n), 2) +
      pow(eff_LooseIsoToLooseIdIso_electronTnP->GetBinContent(n) * err_LooseIdToLooseIdIso_electronTnP->GetBinContent(n),2) +
      pow(eff_LooseIdToLooseIdIso_electronTnP->GetBinContent(n) * err_LooseIsoToLooseIdIso_electronTnP->GetBinContent(n),2)
    ));
    diffErr_Medium_ele->SetBinContent(n, sqrt(
      pow(err_BaselineToMedium_electronTnP->GetBinContent(n), 2) +
      pow(eff_MediumIsoToMediumIdIso_electronTnP->GetBinContent(n) * err_MediumIdToMediumIdIso_electronTnP->GetBinContent(n),2) +
      pow(eff_MediumIdToMediumIdIso_electronTnP->GetBinContent(n) * err_MediumIsoToMediumIdIso_electronTnP->GetBinContent(n),2)
    ));
    diffErr_Tight_ele->SetBinContent(n, sqrt(
      pow(err_BaselineToTight_electronTnP->GetBinContent(n), 2) +
      pow(eff_TightIsoToTightIdIso_electronTnP->GetBinContent(n) * err_TightIdToTightIdIso_electronTnP->GetBinContent(n),2) +
      pow(eff_TightIdToTightIdIso_electronTnP->GetBinContent(n) * err_TightIsoToTightIdIso_electronTnP->GetBinContent(n),2)
    ));
  }}

  // Start drawing stuff 
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;

  // draw |difference| plots
  TCanvas *c_absdiff_Veto_ele = new TCanvas("c_absdiff_Veto_ele","Difference from factorized efficiency, Veto electrons (Data)",800,800);
  absdiff_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Veto_ele->SetTitle("Difference from factorized efficiency, Veto electrons (Data)");
  absdiff_Veto_ele->GetXaxis()->SetTitle("| #eta |");
  absdiff_Veto_ele->GetXaxis()->SetTitleOffset(0.9);
  absdiff_Veto_ele->GetXaxis()->SetTitleSize(0.04);
  absdiff_Veto_ele->GetXaxis()->SetLabelSize(0.02);
  absdiff_Veto_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absdiff_Veto_ele->GetYaxis()->SetTitleOffset(0.9);
  absdiff_Veto_ele->GetYaxis()->SetTitleSize(0.04);
  absdiff_Veto_ele->GetYaxis()->SetLabelSize(0.02);
  absdiff_Veto_ele->GetYaxis()->SetRangeUser(10,100);
  absdiff_Veto_ele->SetMinimum(0);
  absdiff_Veto_ele->SetMaximum(0.2);
  absdiff_Veto_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absdiff_Veto_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absdiff_Veto_ele->Update();
  c_absdiff_Veto_ele->Print((plots_dir+"SingleElectron_absdiff_Veto_ele.png").c_str());
  TCanvas *c_absdiff_Loose_ele = new TCanvas("c_absdiff_Loose_ele","Difference from factorized efficiency, Loose electrons (Data)",800,800);
  absdiff_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Loose_ele->SetTitle("Difference from factorized efficiency, Loose electrons (Data)");
  absdiff_Loose_ele->GetXaxis()->SetTitle("| #eta |");
  absdiff_Loose_ele->GetXaxis()->SetTitleOffset(0.9);
  absdiff_Loose_ele->GetXaxis()->SetTitleSize(0.04);
  absdiff_Loose_ele->GetXaxis()->SetLabelSize(0.02);
  absdiff_Loose_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absdiff_Loose_ele->GetYaxis()->SetTitleOffset(0.9);
  absdiff_Loose_ele->GetYaxis()->SetTitleSize(0.04);
  absdiff_Loose_ele->GetYaxis()->SetLabelSize(0.02);
  absdiff_Loose_ele->GetYaxis()->SetRangeUser(10,100);
  absdiff_Loose_ele->SetMinimum(0);
  absdiff_Loose_ele->SetMaximum(0.2);
  absdiff_Loose_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absdiff_Loose_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absdiff_Loose_ele->Update();
  c_absdiff_Loose_ele->Print((plots_dir+"SingleElectron_absdiff_Loose_ele.png").c_str());
  TCanvas *c_absdiff_Medium_ele = new TCanvas("c_absdiff_Medium_ele","Difference from factorized efficiency, Medium electrons (Data)",800,800);
  absdiff_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Medium_ele->SetTitle("Difference from factorized efficiency, Medium electrons (Data)");
  absdiff_Medium_ele->GetXaxis()->SetTitle("| #eta |");
  absdiff_Medium_ele->GetXaxis()->SetTitleOffset(0.9);
  absdiff_Medium_ele->GetXaxis()->SetTitleSize(0.04);
  absdiff_Medium_ele->GetXaxis()->SetLabelSize(0.02);
  absdiff_Medium_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absdiff_Medium_ele->GetYaxis()->SetTitleOffset(0.9);
  absdiff_Medium_ele->GetYaxis()->SetTitleSize(0.04);
  absdiff_Medium_ele->GetYaxis()->SetLabelSize(0.02);
  absdiff_Medium_ele->GetYaxis()->SetRangeUser(10,100);
  absdiff_Medium_ele->SetMinimum(0);
  absdiff_Medium_ele->SetMaximum(0.2);
  absdiff_Medium_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absdiff_Medium_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absdiff_Medium_ele->Update();
  c_absdiff_Medium_ele->Print((plots_dir+"SingleElectron_absdiff_Medium_ele.png").c_str());
  TCanvas *c_absdiff_Tight_ele = new TCanvas("c_absdiff_Tight_ele","Difference from factorized efficiency, Tight electrons (Data)",800,800);
  absdiff_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Tight_ele->SetTitle("Difference from factorized efficiency, Tight electrons (Data)");
  absdiff_Tight_ele->GetXaxis()->SetTitle("| #eta |");
  absdiff_Tight_ele->GetXaxis()->SetTitleOffset(0.9);
  absdiff_Tight_ele->GetXaxis()->SetTitleSize(0.04);
  absdiff_Tight_ele->GetXaxis()->SetLabelSize(0.02);
  absdiff_Tight_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absdiff_Tight_ele->GetYaxis()->SetTitleOffset(0.9);
  absdiff_Tight_ele->GetYaxis()->SetTitleSize(0.04);
  absdiff_Tight_ele->GetYaxis()->SetLabelSize(0.02);
  absdiff_Tight_ele->GetYaxis()->SetRangeUser(10,100);
  absdiff_Tight_ele->SetMinimum(0);
  absdiff_Tight_ele->SetMaximum(0.2);
  absdiff_Tight_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absdiff_Tight_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absdiff_Tight_ele->Update();
  c_absdiff_Tight_ele->Print((plots_dir+"SingleElectron_absdiff_Tight_ele.png").c_str());

  // Draw unfactorized efficiency plots
  TCanvas *c_eff_Veto_ele = new TCanvas("c_eff_Veto_ele","Unfactorized efficiency, Veto electrons (Data)",800,800);
  eff_BaselineToVeto_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToVeto_electronTnP->SetTitle("Unfactorized efficiency, Veto electrons (Data)");
  eff_BaselineToVeto_electronTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToVeto_electronTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToVeto_electronTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToVeto_electronTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToVeto_electronTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToVeto_electronTnP->SetMinimum(0.0);
  eff_BaselineToVeto_electronTnP->SetMaximum(1);
  eff_BaselineToVeto_electronTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToVeto_electronTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Veto_ele->Update();
  c_eff_Veto_ele->Print((plots_dir+"SingleElectron_eff_Veto_ele.png").c_str());
  TCanvas *c_eff_Loose_ele = new TCanvas("c_eff_Loose_ele","Unfactorized efficiency, Loose electrons (Data)",800,800);
  eff_BaselineToLoose_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToLoose_electronTnP->SetTitle("Unfactorized efficiency, Loose electrons (Data)");
  eff_BaselineToLoose_electronTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToLoose_electronTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToLoose_electronTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToLoose_electronTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToLoose_electronTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToLoose_electronTnP->SetMinimum(0.0);
  eff_BaselineToLoose_electronTnP->SetMaximum(1);
  eff_BaselineToLoose_electronTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToLoose_electronTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Loose_ele->Update();
  c_eff_Loose_ele->Print((plots_dir+"SingleElectron_eff_Loose_ele.png").c_str());
  TCanvas *c_eff_Medium_ele = new TCanvas("c_eff_Medium_ele","Unfactorized efficiency, Medium electrons (Data)",800,800);
  eff_BaselineToMedium_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToMedium_electronTnP->SetTitle("Unfactorized efficiency, Medium electrons (Data)");
  eff_BaselineToMedium_electronTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToMedium_electronTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToMedium_electronTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToMedium_electronTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToMedium_electronTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToMedium_electronTnP->SetMinimum(0.0);
  eff_BaselineToMedium_electronTnP->SetMaximum(1);
  eff_BaselineToMedium_electronTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToMedium_electronTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Medium_ele->Update();
  c_eff_Medium_ele->Print((plots_dir+"SingleElectron_eff_Medium_ele.png").c_str());
  TCanvas *c_eff_Tight_ele = new TCanvas("c_eff_Tight_ele","Unfactorized efficiency, Tight electrons (Data)",800,800);
  eff_BaselineToTight_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToTight_electronTnP->SetTitle("Unfactorized efficiency, Tight electrons (Data)");
  eff_BaselineToTight_electronTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToTight_electronTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToTight_electronTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToTight_electronTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToTight_electronTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToTight_electronTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToTight_electronTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToTight_electronTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToTight_electronTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToTight_electronTnP->SetMinimum(0.0);
  eff_BaselineToTight_electronTnP->SetMaximum(1);
  eff_BaselineToTight_electronTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToTight_electronTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Tight_ele->Update();
  c_eff_Tight_ele->Print((plots_dir+"SingleElectron_eff_Tight_ele.png").c_str());

  // Draw factorized efficiency plots
  TCanvas *c_factEff_Veto_ele = new TCanvas("c_factEff_Veto_ele","Factorized efficiency, Veto electrons (Data)",800,800);
  factEff_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Veto_ele->SetTitle("Factorized efficiency, Veto electrons (Data)");
  factEff_Veto_ele->GetXaxis()->SetTitle("| #eta |");
  factEff_Veto_ele->GetXaxis()->SetTitleOffset(0.9);
  factEff_Veto_ele->GetXaxis()->SetTitleSize(0.04);
  factEff_Veto_ele->GetXaxis()->SetLabelSize(0.02);
  factEff_Veto_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factEff_Veto_ele->GetYaxis()->SetTitleOffset(0.9);
  factEff_Veto_ele->GetYaxis()->SetTitleSize(0.04);
  factEff_Veto_ele->GetYaxis()->SetLabelSize(0.02);
  factEff_Veto_ele->GetYaxis()->SetRangeUser(10,100);
  factEff_Veto_ele->SetMinimum(0.0);
  factEff_Veto_ele->SetMaximum(1);
  factEff_Veto_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factEff_Veto_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factEff_Veto_ele->Update();
  c_factEff_Veto_ele->Print((plots_dir+"SingleElectron_factEff_Veto_ele.png").c_str());
  TCanvas *c_factEff_Loose_ele = new TCanvas("c_factEff_Loose_ele","Factorized efficiency, Loose electrons (Data)",800,800);
  factEff_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Loose_ele->SetTitle("Factorized efficiency, Loose electrons (Data)");
  factEff_Loose_ele->GetXaxis()->SetTitle("| #eta |");
  factEff_Loose_ele->GetXaxis()->SetTitleOffset(0.9);
  factEff_Loose_ele->GetXaxis()->SetTitleSize(0.04);
  factEff_Loose_ele->GetXaxis()->SetLabelSize(0.02);
  factEff_Loose_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factEff_Loose_ele->GetYaxis()->SetTitleOffset(0.9);
  factEff_Loose_ele->GetYaxis()->SetTitleSize(0.04);
  factEff_Loose_ele->GetYaxis()->SetLabelSize(0.02);
  factEff_Loose_ele->GetYaxis()->SetRangeUser(10,100);
  factEff_Loose_ele->SetMinimum(0.0);
  factEff_Loose_ele->SetMaximum(1);
  factEff_Loose_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factEff_Loose_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factEff_Loose_ele->Update();
  c_factEff_Loose_ele->Print((plots_dir+"SingleElectron_factEff_Loose_ele.png").c_str());
  TCanvas *c_factEff_Medium_ele = new TCanvas("c_factEff_Medium_ele","Factorized efficiency, Medium electrons (Data)",800,800);
  factEff_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Medium_ele->SetTitle("Factorized efficiency, Medium electrons (Data)");
  factEff_Medium_ele->GetXaxis()->SetTitle("| #eta |");
  factEff_Medium_ele->GetXaxis()->SetTitleOffset(0.9);
  factEff_Medium_ele->GetXaxis()->SetTitleSize(0.04);
  factEff_Medium_ele->GetXaxis()->SetLabelSize(0.02);
  factEff_Medium_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factEff_Medium_ele->GetYaxis()->SetTitleOffset(0.9);
  factEff_Medium_ele->GetYaxis()->SetTitleSize(0.04);
  factEff_Medium_ele->GetYaxis()->SetLabelSize(0.02);
  factEff_Medium_ele->GetYaxis()->SetRangeUser(10,100);
  factEff_Medium_ele->SetMinimum(0.0);
  factEff_Medium_ele->SetMaximum(1);
  factEff_Medium_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factEff_Medium_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factEff_Medium_ele->Update();
  c_factEff_Medium_ele->Print((plots_dir+"SingleElectron_factEff_Medium_ele.png").c_str());
  TCanvas *c_factEff_Tight_ele = new TCanvas("c_factEff_Tight_ele","Factorized efficiency, Tight electrons (Data)",800,800);
  factEff_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Tight_ele->SetTitle("Factorized efficiency, Tight electrons (Data)");
  factEff_Tight_ele->GetXaxis()->SetTitle("| #eta |");
  factEff_Tight_ele->GetXaxis()->SetTitleOffset(0.9);
  factEff_Tight_ele->GetXaxis()->SetTitleSize(0.04);
  factEff_Tight_ele->GetXaxis()->SetLabelSize(0.02);
  factEff_Tight_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factEff_Tight_ele->GetYaxis()->SetTitleOffset(0.9);
  factEff_Tight_ele->GetYaxis()->SetTitleSize(0.04);
  factEff_Tight_ele->GetYaxis()->SetLabelSize(0.02);
  factEff_Tight_ele->GetYaxis()->SetRangeUser(10,100);
  factEff_Tight_ele->SetMinimum(0.0);
  factEff_Tight_ele->SetMaximum(1);
  factEff_Tight_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factEff_Tight_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factEff_Tight_ele->Update();
  c_factEff_Tight_ele->Print((plots_dir+"SingleElectron_factEff_Tight_ele.png").c_str());
  
  TFile *efficiencies_file = TFile::Open((root_dir+"SingleElectron_efficiencies_electronTnP.root").c_str(),"RECREATE");
  eff_BaselineToVeto_electronTnP   ->Write("eff_Veto_ele");  
  eff_BaselineToLoose_electronTnP  ->Write("eff_Loose_ele"); 
  eff_BaselineToMedium_electronTnP ->Write("eff_Medium_ele"); 
  eff_BaselineToTight_electronTnP  ->Write("eff_Tight_ele");  
  factEff_Veto_ele   ->Write("factEff_Veto_ele");
  factEff_Loose_ele  ->Write("factEff_Loose_ele");
  factEff_Medium_ele ->Write("factEff_Medium_ele");
  factEff_Tight_ele  ->Write("factEff_Tight_ele"); 
  absdiff_Veto_ele   ->Write("absdiff_Veto_ele");
  absdiff_Loose_ele  ->Write("absdiff_Loose_ele");
  absdiff_Medium_ele ->Write("absdiff_Medium_ele");
  absdiff_Tight_ele  ->Write("absdiff_Tight_ele");
  efficiencies_file->Close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void efficiencies_mc_mu(string plots_dir, string root_dir) {
  
  TFile *f_BaselineToLoose_muonTnP          = TFile::Open((plots_dir+"DYJetsToLL_BaselineToLoose_muonTnP/eff.root" ).c_str()       ,"READ");
  TFile *f_BaselineToMedium_muonTnP         = TFile::Open((plots_dir+"DYJetsToLL_BaselineToMedium_muonTnP/eff.root").c_str()       ,"READ");
  TFile *f_BaselineToTight_muonTnP          = TFile::Open((plots_dir+"DYJetsToLL_BaselineToTight_muonTnP/eff.root" ).c_str()       ,"READ");
  TFile *f_BaselineToVeto_muonTnP           = TFile::Open((plots_dir+"DYJetsToLL_BaselineToVeto_muonTnP/eff.root"  ).c_str()       ,"READ");
  
  TH2D * eff_BaselineToLoose_muonTnP        = (TH2D*) f_BaselineToLoose_muonTnP          ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToMedium_muonTnP       = (TH2D*) f_BaselineToMedium_muonTnP         ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToTight_muonTnP        = (TH2D*) f_BaselineToTight_muonTnP          ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToVeto_muonTnP         = (TH2D*) f_BaselineToVeto_muonTnP           ->Get("hEffEtaPt"); 
  
  TH2D * err_BaselineToLoose_muonTnP        = (TH2D*) f_BaselineToLoose_muonTnP          ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToMedium_muonTnP       = (TH2D*) f_BaselineToMedium_muonTnP         ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToTight_muonTnP        = (TH2D*) f_BaselineToTight_muonTnP          ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToVeto_muonTnP         = (TH2D*) f_BaselineToVeto_muonTnP           ->Get("hErrlEtaPt"); 

  for(int i_eta = 1; i_eta <= n_mu_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_mu_pt_bins; i_pt++) {
    int n = eff_BaselineToLoose_muonTnP->GetBin(i_eta, i_pt);
    eff_BaselineToLoose_muonTnP        ->SetBinError(n, err_BaselineToLoose_muonTnP        ->GetBinContent(n) ); 
    eff_BaselineToMedium_muonTnP       ->SetBinError(n, err_BaselineToMedium_muonTnP       ->GetBinContent(n) );
    eff_BaselineToTight_muonTnP        ->SetBinError(n, err_BaselineToTight_muonTnP        ->GetBinContent(n) );
    eff_BaselineToVeto_muonTnP         ->SetBinError(n, err_BaselineToVeto_muonTnP         ->GetBinContent(n) );
  }}
  
  // Start drawing stuff 
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;

  // Draw unfactorized efficiency plots
  TCanvas *c_eff_Veto_mu = new TCanvas("c_eff_Veto_mu","Unfactorized Veto muon efficiency (MC)",800,800);
  eff_BaselineToVeto_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToVeto_muonTnP->SetTitle("Unfactorized Veto muon efficiency (MC)");
  eff_BaselineToVeto_muonTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToVeto_muonTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToVeto_muonTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToVeto_muonTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToVeto_muonTnP->SetMinimum(0.4);
  eff_BaselineToVeto_muonTnP->SetMaximum(1);
  eff_BaselineToVeto_muonTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToVeto_muonTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Veto_mu->Update();
  c_eff_Veto_mu->Print((plots_dir+"DYJetsToLL_eff_Veto_mu.png").c_str());
  TCanvas *c_eff_Loose_mu = new TCanvas("c_eff_Loose_mu","Unfactorized Loose muon efficiency (MC)",800,800);
  eff_BaselineToLoose_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToLoose_muonTnP->SetTitle("Unfactorized Loose muon efficiency (MC)");
  eff_BaselineToLoose_muonTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToLoose_muonTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToLoose_muonTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToLoose_muonTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToLoose_muonTnP->SetMinimum(0.4);
  eff_BaselineToLoose_muonTnP->SetMaximum(1);
  eff_BaselineToLoose_muonTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToLoose_muonTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Loose_mu->Update();
  c_eff_Loose_mu->Print((plots_dir+"DYJetsToLL_eff_Loose_mu.png").c_str());
  TCanvas *c_eff_Medium_mu = new TCanvas("c_eff_Medium_mu","Unfactorized Medium muon efficiency (MC)",800,800);
  eff_BaselineToMedium_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToMedium_muonTnP->SetTitle("Unfactorized Medium muon efficiency (MC)");
  eff_BaselineToMedium_muonTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToMedium_muonTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToMedium_muonTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToMedium_muonTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToMedium_muonTnP->SetMinimum(0.4);
  eff_BaselineToMedium_muonTnP->SetMaximum(1);
  eff_BaselineToMedium_muonTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToMedium_muonTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Medium_mu->Update();
  c_eff_Medium_mu->Print((plots_dir+"DYJetsToLL_eff_Medium_mu.png").c_str());
  TCanvas *c_eff_Tight_mu = new TCanvas("c_eff_Tight_mu","Unfactorized Tight muon efficiency (MC)",800,800);
  eff_BaselineToTight_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToTight_muonTnP->SetTitle("Unfactorized Tight muon efficiency (MC)");
  eff_BaselineToTight_muonTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToTight_muonTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToTight_muonTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToTight_muonTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToTight_muonTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToTight_muonTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToTight_muonTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToTight_muonTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToTight_muonTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToTight_muonTnP->SetMinimum(0.4);
  eff_BaselineToTight_muonTnP->SetMaximum(1);
  eff_BaselineToTight_muonTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToTight_muonTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Tight_mu->Update();
  c_eff_Tight_mu->Print((plots_dir+"DYJetsToLL_eff_Tight_mu.png").c_str());

  TFile *efficiencies_file = TFile::Open((root_dir+"DYJetsToLL_efficiencies_muonTnP.root").c_str() ,"RECREATE");
  eff_BaselineToVeto_muonTnP   ->Write("eff_Veto_mu");  
  eff_BaselineToLoose_muonTnP  ->Write("eff_Loose_mu"); 
  eff_BaselineToMedium_muonTnP ->Write("eff_Medium_mu"); 
  eff_BaselineToTight_muonTnP  ->Write("eff_Tight_mu");  
  efficiencies_file->Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void efficiencies_data_mu(string plots_dir, string root_dir) {
  
  TFile *f_BaselineToLoose_muonTnP          = TFile::Open((plots_dir+"SingleMuon_BaselineToLoose_muonTnP/eff.root" ).c_str()       ,"READ");
  TFile *f_BaselineToMedium_muonTnP         = TFile::Open((plots_dir+"SingleMuon_BaselineToMedium_muonTnP/eff.root").c_str()       ,"READ");
  TFile *f_BaselineToTight_muonTnP          = TFile::Open((plots_dir+"SingleMuon_BaselineToTight_muonTnP/eff.root" ).c_str()       ,"READ");
  TFile *f_BaselineToVeto_muonTnP           = TFile::Open((plots_dir+"SingleMuon_BaselineToVeto_muonTnP/eff.root"  ).c_str()       ,"READ");
  
  TH2D * eff_BaselineToLoose_muonTnP        = (TH2D*) f_BaselineToLoose_muonTnP          ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToMedium_muonTnP       = (TH2D*) f_BaselineToMedium_muonTnP         ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToTight_muonTnP        = (TH2D*) f_BaselineToTight_muonTnP          ->Get("hEffEtaPt"); 
  TH2D * eff_BaselineToVeto_muonTnP         = (TH2D*) f_BaselineToVeto_muonTnP           ->Get("hEffEtaPt"); 
  
  TH2D * err_BaselineToLoose_muonTnP        = (TH2D*) f_BaselineToLoose_muonTnP          ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToMedium_muonTnP       = (TH2D*) f_BaselineToMedium_muonTnP         ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToTight_muonTnP        = (TH2D*) f_BaselineToTight_muonTnP          ->Get("hErrlEtaPt"); 
  TH2D * err_BaselineToVeto_muonTnP         = (TH2D*) f_BaselineToVeto_muonTnP           ->Get("hErrlEtaPt"); 

  for(int i_eta = 1; i_eta <= n_mu_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_mu_pt_bins; i_pt++) {
    int n = eff_BaselineToLoose_muonTnP->GetBin(i_eta, i_pt);
    eff_BaselineToLoose_muonTnP        ->SetBinError(n, err_BaselineToLoose_muonTnP        ->GetBinContent(n) ); 
    eff_BaselineToMedium_muonTnP       ->SetBinError(n, err_BaselineToMedium_muonTnP       ->GetBinContent(n) );
    eff_BaselineToTight_muonTnP        ->SetBinError(n, err_BaselineToTight_muonTnP        ->GetBinContent(n) );
    eff_BaselineToVeto_muonTnP         ->SetBinError(n, err_BaselineToVeto_muonTnP         ->GetBinContent(n) );
  }}
  
  // Start drawing stuff 
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;

  // Draw unfactorized efficiency plots
  TCanvas *c_eff_Veto_mu = new TCanvas("c_eff_Veto_mu","Unfactorized Veto muon efficiency (Data)",800,800);
  eff_BaselineToVeto_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToVeto_muonTnP->SetTitle("Unfactorized Veto muon efficiency (Data)");
  eff_BaselineToVeto_muonTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToVeto_muonTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToVeto_muonTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToVeto_muonTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToVeto_muonTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToVeto_muonTnP->SetMinimum(0.4);
  eff_BaselineToVeto_muonTnP->SetMaximum(1);
  eff_BaselineToVeto_muonTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToVeto_muonTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Veto_mu->Update();
  c_eff_Veto_mu->Print((plots_dir+"SingleMuon_eff_Veto_mu.png").c_str());
  TCanvas *c_eff_Loose_mu = new TCanvas("c_eff_Loose_mu","Unfactorized Loose muon efficiency (Data)",800,800);
  eff_BaselineToLoose_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToLoose_muonTnP->SetTitle("Unfactorized Loose muon efficiency (Data)");
  eff_BaselineToLoose_muonTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToLoose_muonTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToLoose_muonTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToLoose_muonTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToLoose_muonTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToLoose_muonTnP->SetMinimum(0.4);
  eff_BaselineToLoose_muonTnP->SetMaximum(1);
  eff_BaselineToLoose_muonTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToLoose_muonTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Loose_mu->Update();
  c_eff_Loose_mu->Print((plots_dir+"SingleMuon_eff_Loose_mu.png").c_str());
  TCanvas *c_eff_Medium_mu = new TCanvas("c_eff_Medium_mu","Unfactorized Medium muon efficiency (Data)",800,800);
  eff_BaselineToMedium_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToMedium_muonTnP->SetTitle("Unfactorized Medium muon efficiency (Data)");
  eff_BaselineToMedium_muonTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToMedium_muonTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToMedium_muonTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToMedium_muonTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToMedium_muonTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToMedium_muonTnP->SetMinimum(0.4);
  eff_BaselineToMedium_muonTnP->SetMaximum(1);
  eff_BaselineToMedium_muonTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToMedium_muonTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Medium_mu->Update();
  c_eff_Medium_mu->Print((plots_dir+"SingleMuon_eff_Medium_mu.png").c_str());
  TCanvas *c_eff_Tight_mu = new TCanvas("c_eff_Tight_mu","Unfactorized Tight muon efficiency (Data)",800,800);
  eff_BaselineToTight_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToTight_muonTnP->SetTitle("Unfactorized Tight muon efficiency (Data)");
  eff_BaselineToTight_muonTnP->GetXaxis()->SetTitle("| #eta |");
  eff_BaselineToTight_muonTnP->GetXaxis()->SetTitleOffset(0.9);
  eff_BaselineToTight_muonTnP->GetXaxis()->SetTitleSize(0.04);
  eff_BaselineToTight_muonTnP->GetXaxis()->SetLabelSize(0.02);
  eff_BaselineToTight_muonTnP->GetYaxis()->SetTitle("p_{T} [GeV]");
  eff_BaselineToTight_muonTnP->GetYaxis()->SetTitleOffset(0.9);
  eff_BaselineToTight_muonTnP->GetYaxis()->SetTitleSize(0.04);
  eff_BaselineToTight_muonTnP->GetYaxis()->SetLabelSize(0.02);
  eff_BaselineToTight_muonTnP->GetYaxis()->SetRangeUser(10,100);
  eff_BaselineToTight_muonTnP->SetMinimum(0.4);
  eff_BaselineToTight_muonTnP->SetMaximum(1);
  eff_BaselineToTight_muonTnP->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) eff_BaselineToTight_muonTnP->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_eff_Tight_mu->Update();
  c_eff_Tight_mu->Print((plots_dir+"SingleMuon_eff_Tight_mu.png").c_str());

  TFile *efficiencies_file = TFile::Open((root_dir+"SingleMuon_efficiencies_muonTnP.root").c_str(),"RECREATE");
  eff_BaselineToVeto_muonTnP   ->Write("eff_Veto_mu");  
  eff_BaselineToLoose_muonTnP  ->Write("eff_Loose_mu"); 
  eff_BaselineToMedium_muonTnP ->Write("eff_Medium_mu"); 
  eff_BaselineToTight_muonTnP  ->Write("eff_Tight_mu");  
  
  efficiencies_file->Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void make_ele_scale_factors(string plots_dir, string root_dir) {

  TFile *f_mc   = TFile::Open((root_dir+"DYJetsToLL_efficiencies_electronTnP.root").c_str(),"READ");
  TFile *f_data = TFile::Open((root_dir+"SingleElectron_efficiencies_electronTnP.root").c_str(),"READ");

  // compute Veto scalefactors
  TH2D *eff_Veto_ele_data   = (TH2D*) f_data ->Get("eff_Veto_ele");
  TH2D *eff_Veto_ele_mc     = (TH2D*) f_mc   ->Get("eff_Veto_ele");
  TH2D *unfactorized_scalefactors_Veto_ele = new TH2D(
    "unfactorized_scalefactors_Veto_ele",
    "Electron Veto ID scale factors from Unfactorized efficiencies",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  (*unfactorized_scalefactors_Veto_ele) = (*eff_Veto_ele_data) / (*eff_Veto_ele_mc);
  TH2D *factEff_Veto_ele_data   = (TH2D*) f_data ->Get("factEff_Veto_ele");
  TH2D *factEff_Veto_ele_mc     = (TH2D*) f_mc   ->Get("factEff_Veto_ele");
  TH2D *factorized_scalefactors_Veto_ele = new TH2D(
    "factorized_scalefactors_Veto_ele",
    "Electron Veto ID scale factors from Factorized efficiencies",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  (*factorized_scalefactors_Veto_ele) = (*factEff_Veto_ele_data) / (*factEff_Veto_ele_mc);

  // compute Loose scalefactors
  TH2D *eff_Loose_ele_data   = (TH2D*) f_data ->Get("eff_Loose_ele");
  TH2D *eff_Loose_ele_mc     = (TH2D*) f_mc   ->Get("eff_Loose_ele");
  TH2D *unfactorized_scalefactors_Loose_ele = new TH2D(
    "unfactorized_scalefactors_Loose_ele",
    "Electron Loose ID scale factors from Unfactorized efficiencies",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  (*unfactorized_scalefactors_Loose_ele) = (*eff_Loose_ele_data) / (*eff_Loose_ele_mc);
  TH2D *factEff_Loose_ele_data   = (TH2D*) f_data ->Get("factEff_Loose_ele");
  TH2D *factEff_Loose_ele_mc     = (TH2D*) f_mc   ->Get("factEff_Loose_ele");
  TH2D *factorized_scalefactors_Loose_ele = new TH2D(
    "factorized_scalefactors_Loose_ele",
    "Electron Loose ID scale factors from Factorized efficiencies",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  (*factorized_scalefactors_Loose_ele) = (*factEff_Loose_ele_data) / (*factEff_Loose_ele_mc);

  // compute Medium scalefactors
  TH2D *eff_Medium_ele_data   = (TH2D*) f_data ->Get("eff_Medium_ele");
  TH2D *eff_Medium_ele_mc     = (TH2D*) f_mc   ->Get("eff_Medium_ele");
  TH2D *unfactorized_scalefactors_Medium_ele = new TH2D(
    "unfactorized_scalefactors_Medium_ele",
    "Electron Medium ID scale factors from Unfactorized efficiencies",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  (*unfactorized_scalefactors_Medium_ele) = (*eff_Medium_ele_data) / (*eff_Medium_ele_mc);
  TH2D *factEff_Medium_ele_data   = (TH2D*) f_data ->Get("factEff_Medium_ele");
  TH2D *factEff_Medium_ele_mc     = (TH2D*) f_mc   ->Get("factEff_Medium_ele");
  TH2D *factorized_scalefactors_Medium_ele = new TH2D(
    "factorized_scalefactors_Medium_ele",
    "Electron Medium ID scale factors from Factorized efficiencies",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  (*factorized_scalefactors_Medium_ele) = (*factEff_Medium_ele_data) / (*factEff_Medium_ele_mc);

  // compute Tight scalefactors
  TH2D *eff_Tight_ele_data   = (TH2D*) f_data ->Get("eff_Tight_ele");
  TH2D *eff_Tight_ele_mc     = (TH2D*) f_mc   ->Get("eff_Tight_ele");
  TH2D *unfactorized_scalefactors_Tight_ele = new TH2D(
    "unfactorized_scalefactors_Tight_ele",
    "Electron Tight ID scale factors from Unfactorized efficiencies",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  (*unfactorized_scalefactors_Tight_ele) = (*eff_Tight_ele_data) / (*eff_Tight_ele_mc);
  TH2D *factEff_Tight_ele_data   = (TH2D*) f_data ->Get("factEff_Tight_ele");
  TH2D *factEff_Tight_ele_mc     = (TH2D*) f_mc   ->Get("factEff_Tight_ele");
  TH2D *factorized_scalefactors_Tight_ele = new TH2D(
    "factorized_scalefactors_Tight_ele",
    "Electron Tight ID scale factors from Factorized efficiencies",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  (*factorized_scalefactors_Tight_ele) = (*factEff_Tight_ele_data) / (*factEff_Tight_ele_mc);
  
  // Start drawing stuff 
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;
  TLine *oneline=new TLine(10,1,100,1);
  oneline->SetLineColor(1);
  oneline->SetLineStyle(3);
  Int_t mit_red  = 1861; 
  Int_t mit_gray = 1862; 
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);

  // Draw unfactorized scalefactor plots
  TCanvas *c_unfactorized_scalefactors_Veto_ele = new TCanvas("c_unfactorized_scalefactors_Veto_ele","Unfactorized Veto electron scale factors (Data/MC)",800,800);
  unfactorized_scalefactors_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Veto_ele->SetTitle("Unfactorized Veto electron scale factors (Data/MC)");
  unfactorized_scalefactors_Veto_ele->GetXaxis()->SetTitle("| #eta |");
  unfactorized_scalefactors_Veto_ele->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Veto_ele->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Veto_ele->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Veto_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Veto_ele->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Veto_ele->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Veto_ele->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Veto_ele->GetYaxis()->SetRangeUser(10,100);
  unfactorized_scalefactors_Veto_ele->SetMinimum(0.7);
  unfactorized_scalefactors_Veto_ele->SetMaximum(1.3);
  unfactorized_scalefactors_Veto_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) unfactorized_scalefactors_Veto_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_unfactorized_scalefactors_Veto_ele->Update();
  c_unfactorized_scalefactors_Veto_ele->Print((plots_dir+"unfactorized_scalefactors_Veto_ele.png").c_str());
  TCanvas *c_unfactorized_scalefactors_Loose_ele = new TCanvas("c_unfactorized_scalefactors_Loose_ele","Unfactorized Loose electron scale factors (Data/MC)",800,800);
  unfactorized_scalefactors_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Loose_ele->SetTitle("Unfactorized Medium electron scale factors, Electron Loose ID (Data/MC)");
  unfactorized_scalefactors_Loose_ele->GetXaxis()->SetTitle("| #eta |");
  unfactorized_scalefactors_Loose_ele->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Loose_ele->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Loose_ele->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Loose_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Loose_ele->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Loose_ele->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Loose_ele->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Loose_ele->GetYaxis()->SetRangeUser(10,100);
  unfactorized_scalefactors_Loose_ele->SetMinimum(0.7);
  unfactorized_scalefactors_Loose_ele->SetMaximum(1.3);
  unfactorized_scalefactors_Loose_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) unfactorized_scalefactors_Loose_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_unfactorized_scalefactors_Loose_ele->Update();
  c_unfactorized_scalefactors_Loose_ele->Print((plots_dir+"unfactorized_scalefactors_Loose_ele.png").c_str());
  TCanvas *c_unfactorized_scalefactors_Medium_ele = new TCanvas("c_unfactorized_scalefactors_Medium_ele","Unfactorized Medium electron scale factors (Data/MC)",800,800);
  unfactorized_scalefactors_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Medium_ele->SetTitle("Unfactorized Medium electron scale factors (Data/MC)");
  unfactorized_scalefactors_Medium_ele->GetXaxis()->SetTitle("| #eta |");
  unfactorized_scalefactors_Medium_ele->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Medium_ele->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Medium_ele->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Medium_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Medium_ele->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Medium_ele->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Medium_ele->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Medium_ele->GetYaxis()->SetRangeUser(10,100);
  unfactorized_scalefactors_Medium_ele->SetMinimum(0.7);
  unfactorized_scalefactors_Medium_ele->SetMaximum(1.3);
  unfactorized_scalefactors_Medium_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) unfactorized_scalefactors_Medium_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_unfactorized_scalefactors_Medium_ele->Update();
  c_unfactorized_scalefactors_Medium_ele->Print((plots_dir+"unfactorized_scalefactors_Medium_ele.png").c_str());
  TCanvas *c_unfactorized_scalefactors_Tight_ele = new TCanvas("c_unfactorized_scalefactors_Tight_ele","Unfactorized Tight electron scale factors (Data/MC)",800,800);
  unfactorized_scalefactors_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Tight_ele->SetTitle("Unfactorized Tight elecron scale factors (Data/MC)");
  unfactorized_scalefactors_Tight_ele->GetXaxis()->SetTitle("| #eta |");
  unfactorized_scalefactors_Tight_ele->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Tight_ele->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Tight_ele->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Tight_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Tight_ele->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Tight_ele->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Tight_ele->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Tight_ele->GetYaxis()->SetRangeUser(10,100);
  unfactorized_scalefactors_Tight_ele->SetMinimum(0.7);
  unfactorized_scalefactors_Tight_ele->SetMaximum(1.3);
  unfactorized_scalefactors_Tight_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) unfactorized_scalefactors_Tight_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_unfactorized_scalefactors_Tight_ele->Update();
  c_unfactorized_scalefactors_Tight_ele->Print((plots_dir+"unfactorized_scalefactors_Tight_ele.png").c_str());

  // Draw factorized scalefactor plots
  TCanvas *c_factorized_scalefactors_Veto_ele = new TCanvas("c_factorized_scalefactors_Veto_ele","Factorized Veto electron scale factors (Data/MC)",800,800);
  factorized_scalefactors_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factorized_scalefactors_Veto_ele->SetTitle("Factorized Veto electron scale factors (Data/MC)");
  factorized_scalefactors_Veto_ele->GetXaxis()->SetTitle("| #eta |");
  factorized_scalefactors_Veto_ele->GetXaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Veto_ele->GetXaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Veto_ele->GetXaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Veto_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factorized_scalefactors_Veto_ele->GetYaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Veto_ele->GetYaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Veto_ele->GetYaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Veto_ele->GetYaxis()->SetRangeUser(10,100);
  factorized_scalefactors_Veto_ele->SetMinimum(0.7);
  factorized_scalefactors_Veto_ele->SetMaximum(1.3);
  factorized_scalefactors_Veto_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factorized_scalefactors_Veto_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factorized_scalefactors_Veto_ele->Update();
  c_factorized_scalefactors_Veto_ele->Print((plots_dir+"factorized_scalefactors_Veto_ele.png").c_str());
  TCanvas *c_factorized_scalefactors_Loose_ele = new TCanvas("c_factorized_scalefactors_Loose_ele","Factorized Loose electron scale factors (Data/MC)",800,800);
  factorized_scalefactors_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factorized_scalefactors_Loose_ele->SetTitle("Factorized Loose electron scale factors (Data/MC)");
  factorized_scalefactors_Loose_ele->GetXaxis()->SetTitle("| #eta |");
  factorized_scalefactors_Loose_ele->GetXaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Loose_ele->GetXaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Loose_ele->GetXaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Loose_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factorized_scalefactors_Loose_ele->GetYaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Loose_ele->GetYaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Loose_ele->GetYaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Loose_ele->GetYaxis()->SetRangeUser(10,100);
  factorized_scalefactors_Loose_ele->SetMinimum(0.7);
  factorized_scalefactors_Loose_ele->SetMaximum(1.3);
  factorized_scalefactors_Loose_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factorized_scalefactors_Loose_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factorized_scalefactors_Loose_ele->Update();
  c_factorized_scalefactors_Loose_ele->Print((plots_dir+"factorized_scalefactors_Loose_ele.png").c_str());
  TCanvas *c_factorized_scalefactors_Medium_ele = new TCanvas("c_factorized_scalefactors_Medium_ele","Factorized Medium electron scale factors (Data/MC)",800,800);
  factorized_scalefactors_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factorized_scalefactors_Medium_ele->SetTitle("Factorized Medium electron scale factors (Data/MC)");
  factorized_scalefactors_Medium_ele->GetXaxis()->SetTitle("| #eta |");
  factorized_scalefactors_Medium_ele->GetXaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Medium_ele->GetXaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Medium_ele->GetXaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Medium_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factorized_scalefactors_Medium_ele->GetYaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Medium_ele->GetYaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Medium_ele->GetYaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Medium_ele->GetYaxis()->SetRangeUser(10,100);
  factorized_scalefactors_Medium_ele->SetMinimum(0.7);
  factorized_scalefactors_Medium_ele->SetMaximum(1.3);
  factorized_scalefactors_Medium_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factorized_scalefactors_Medium_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factorized_scalefactors_Medium_ele->Update();
  c_factorized_scalefactors_Medium_ele->Print((plots_dir+"factorized_scalefactors_Medium_ele.png").c_str());
  TCanvas *c_factorized_scalefactors_Tight_ele = new TCanvas("c_factorized_scalefactors_Tight_ele","Factorized Tight electron scale factors (Data/MC)",800,800);
  factorized_scalefactors_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factorized_scalefactors_Tight_ele->SetTitle("Factorized Tight electron scale factors (Data/MC)");
  factorized_scalefactors_Tight_ele->GetXaxis()->SetTitle("| #eta |");
  factorized_scalefactors_Tight_ele->GetXaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Tight_ele->GetXaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Tight_ele->GetXaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Tight_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  factorized_scalefactors_Tight_ele->GetYaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Tight_ele->GetYaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Tight_ele->GetYaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Tight_ele->GetYaxis()->SetRangeUser(10,100);
  factorized_scalefactors_Tight_ele->SetMinimum(0.7);
  factorized_scalefactors_Tight_ele->SetMaximum(1.3);
  factorized_scalefactors_Tight_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) factorized_scalefactors_Tight_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_factorized_scalefactors_Tight_ele->Update();
  c_factorized_scalefactors_Tight_ele->Print((plots_dir+"factorized_scalefactors_Tight_ele.png").c_str());

  TCanvas *c_factorized_scalefactors_Veto_ele_1d = new TCanvas("c_factorized_scalefactors_Veto_ele_1d","Factorized Veto electron scale factors",800,600);
  TH1D *factorized_scalefactors_Veto_ele_eta1 = factorized_scalefactors_Veto_ele->ProjectionY("factorized_scalefactors_Veto_ele_eta1", 1, 1);
  TH1D *factorized_scalefactors_Veto_ele_eta2 = factorized_scalefactors_Veto_ele->ProjectionY("factorized_scalefactors_Veto_ele_eta2", 2, 2);
  factorized_scalefactors_Veto_ele_eta1->Draw("P0 E1");
  factorized_scalefactors_Veto_ele_eta1->SetTitle("Factorized Veto electron scale factors");
  factorized_scalefactors_Veto_ele_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
  factorized_scalefactors_Veto_ele_eta1->GetXaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Veto_ele_eta1->GetXaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Veto_ele_eta1->GetXaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Veto_ele_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
  factorized_scalefactors_Veto_ele_eta1->GetYaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Veto_ele_eta1->GetYaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Veto_ele_eta1->GetYaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Veto_ele_eta1->SetMinimum(0.8);
  factorized_scalefactors_Veto_ele_eta1->SetMaximum(1.2);
  factorized_scalefactors_Veto_ele_eta1->SetMarkerStyle(20);
  factorized_scalefactors_Veto_ele_eta1->SetMarkerColor(1);
  factorized_scalefactors_Veto_ele_eta1->SetLineColor(1);
  factorized_scalefactors_Veto_ele_eta2->SetMarkerStyle(21);
  factorized_scalefactors_Veto_ele_eta2->SetMarkerColor(1861);
  factorized_scalefactors_Veto_ele_eta2->SetLineColor(1861);
  factorized_scalefactors_Veto_ele_eta2->Draw("P0 E1 SAME");
  oneline->Draw("SAME");
  TLegend *legend_Veto_1d=new TLegend(.7,.7,.85,.85);
  legend_Veto_1d->AddEntry(factorized_scalefactors_Veto_ele_eta1,"| #eta | < 1.479", "lp");
  legend_Veto_1d->AddEntry(factorized_scalefactors_Veto_ele_eta2,"| #eta | > 1.479", "lp");
  legend_Veto_1d->SetFillColor(0);
  legend_Veto_1d->Draw("SAME");
  c_factorized_scalefactors_Veto_ele_1d->Update();
  c_factorized_scalefactors_Veto_ele_1d->Print((plots_dir+"factorized_scalefactors_Veto_ele_1d.png").c_str());
  
  TCanvas *c_factorized_scalefactors_Loose_ele_1d = new TCanvas("c_factorized_scalefactors_Loose_ele_1d","Factorized Loose electron scale factors",800,600);
  TH1D *factorized_scalefactors_Loose_ele_eta1 = factorized_scalefactors_Loose_ele->ProjectionY("factorized_scalefactors_Loose_ele_eta1", 1, 1);
  TH1D *factorized_scalefactors_Loose_ele_eta2 = factorized_scalefactors_Loose_ele->ProjectionY("factorized_scalefactors_Loose_ele_eta2", 2, 2);
  factorized_scalefactors_Loose_ele_eta1->Draw("P0 E1");
  factorized_scalefactors_Loose_ele_eta1->SetTitle("Factorized Loose electron scale factors");
  factorized_scalefactors_Loose_ele_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
  factorized_scalefactors_Loose_ele_eta1->GetXaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Loose_ele_eta1->GetXaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Loose_ele_eta1->GetXaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Loose_ele_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
  factorized_scalefactors_Loose_ele_eta1->GetYaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Loose_ele_eta1->GetYaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Loose_ele_eta1->GetYaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Loose_ele_eta1->SetMinimum(0.8);
  factorized_scalefactors_Loose_ele_eta1->SetMaximum(1.2);
  factorized_scalefactors_Loose_ele_eta1->SetMarkerStyle(20);
  factorized_scalefactors_Loose_ele_eta1->SetMarkerColor(1);
  factorized_scalefactors_Loose_ele_eta1->SetLineColor(1);
  factorized_scalefactors_Loose_ele_eta2->SetMarkerStyle(21);
  factorized_scalefactors_Loose_ele_eta2->SetMarkerColor(1861);
  factorized_scalefactors_Loose_ele_eta2->SetLineColor(1861);
  factorized_scalefactors_Loose_ele_eta2->Draw("P0 E1 SAME");
  oneline->Draw("SAME");
  TLegend *legend_Loose_1d=new TLegend(.7,.7,.85,.85);
  legend_Loose_1d->AddEntry(factorized_scalefactors_Loose_ele_eta1,"| #eta | < 1.479", "lp");
  legend_Loose_1d->AddEntry(factorized_scalefactors_Loose_ele_eta2,"| #eta | > 1.479", "lp");
  legend_Loose_1d->SetFillColor(0);
  legend_Loose_1d->Draw("SAME");
  c_factorized_scalefactors_Loose_ele_1d->Update();
  c_factorized_scalefactors_Loose_ele_1d->Print((plots_dir+"factorized_scalefactors_Loose_ele_1d.png").c_str());
  
  TCanvas *c_factorized_scalefactors_Medium_ele_1d = new TCanvas("c_factorized_scalefactors_Medium_ele_1d","Factorized Medium electron scale factors",800,600);
  TH1D *factorized_scalefactors_Medium_ele_eta1 = factorized_scalefactors_Medium_ele->ProjectionY("factorized_scalefactors_Medium_ele_eta1", 1, 1);
  TH1D *factorized_scalefactors_Medium_ele_eta2 = factorized_scalefactors_Medium_ele->ProjectionY("factorized_scalefactors_Medium_ele_eta2", 2, 2);
  factorized_scalefactors_Medium_ele_eta1->Draw("P0 E1");
  factorized_scalefactors_Medium_ele_eta1->SetTitle("Factorized Medium electron scale factors");
  factorized_scalefactors_Medium_ele_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
  factorized_scalefactors_Medium_ele_eta1->GetXaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Medium_ele_eta1->GetXaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Medium_ele_eta1->GetXaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Medium_ele_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
  factorized_scalefactors_Medium_ele_eta1->GetYaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Medium_ele_eta1->GetYaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Medium_ele_eta1->GetYaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Medium_ele_eta1->SetMinimum(0.8);
  factorized_scalefactors_Medium_ele_eta1->SetMaximum(1.2);
  factorized_scalefactors_Medium_ele_eta1->SetMarkerStyle(20);
  factorized_scalefactors_Medium_ele_eta1->SetMarkerColor(1);
  factorized_scalefactors_Medium_ele_eta1->SetLineColor(1);
  factorized_scalefactors_Medium_ele_eta2->SetMarkerStyle(21);
  factorized_scalefactors_Medium_ele_eta2->SetMarkerColor(1861);
  factorized_scalefactors_Medium_ele_eta2->SetLineColor(1861);
  factorized_scalefactors_Medium_ele_eta2->Draw("P0 E1 SAME");
  oneline->Draw("SAME");
  TLegend *legend_Medium_1d=new TLegend(.7,.7,.85,.85);
  legend_Medium_1d->AddEntry(factorized_scalefactors_Medium_ele_eta1,"| #eta | < 1.479", "lp");
  legend_Medium_1d->AddEntry(factorized_scalefactors_Medium_ele_eta2,"| #eta | > 1.479", "lp");
  legend_Medium_1d->SetFillColor(0);
  legend_Medium_1d->Draw("SAME");
  c_factorized_scalefactors_Medium_ele_1d->Update();
  c_factorized_scalefactors_Medium_ele_1d->Print((plots_dir+"factorized_scalefactors_Medium_ele_1d.png").c_str());
  
  TCanvas *c_factorized_scalefactors_Tight_ele_1d = new TCanvas("c_factorized_scalefactors_Tight_ele_1d","Factorized Tight electron scale factors",800,600);
  TH1D *factorized_scalefactors_Tight_ele_eta1 = factorized_scalefactors_Tight_ele->ProjectionY("factorized_scalefactors_Tight_ele_eta1", 1, 1);
  TH1D *factorized_scalefactors_Tight_ele_eta2 = factorized_scalefactors_Tight_ele->ProjectionY("factorized_scalefactors_Tight_ele_eta2", 2, 2);
  factorized_scalefactors_Tight_ele_eta1->Draw("P0 E1");
  factorized_scalefactors_Tight_ele_eta1->SetTitle("Factorized Tight electron scale factors");
  factorized_scalefactors_Tight_ele_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
  factorized_scalefactors_Tight_ele_eta1->GetXaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Tight_ele_eta1->GetXaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Tight_ele_eta1->GetXaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Tight_ele_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
  factorized_scalefactors_Tight_ele_eta1->GetYaxis()->SetTitleOffset(0.9);
  factorized_scalefactors_Tight_ele_eta1->GetYaxis()->SetTitleSize(0.04);
  factorized_scalefactors_Tight_ele_eta1->GetYaxis()->SetLabelSize(0.02);
  factorized_scalefactors_Tight_ele_eta1->SetMinimum(0.8);
  factorized_scalefactors_Tight_ele_eta1->SetMaximum(1.2);
  factorized_scalefactors_Tight_ele_eta1->SetMarkerStyle(20);
  factorized_scalefactors_Tight_ele_eta1->SetMarkerColor(1);
  factorized_scalefactors_Tight_ele_eta1->SetLineColor(1);
  factorized_scalefactors_Tight_ele_eta2->SetMarkerStyle(21);
  factorized_scalefactors_Tight_ele_eta2->SetMarkerColor(1861);
  factorized_scalefactors_Tight_ele_eta2->SetLineColor(1861);
  factorized_scalefactors_Tight_ele_eta2->Draw("P0 E1 SAME");
  oneline->Draw("SAME");
  TLegend *legend_Tight_1d=new TLegend(.7,.7,.85,.85);
  legend_Tight_1d->AddEntry(factorized_scalefactors_Tight_ele_eta1,"| #eta | < 1.479", "lp");
  legend_Tight_1d->AddEntry(factorized_scalefactors_Tight_ele_eta2,"| #eta | > 1.479", "lp");
  legend_Tight_1d->SetFillColor(0);
  legend_Tight_1d->Draw("SAME");
  c_factorized_scalefactors_Tight_ele_1d->Update();
  c_factorized_scalefactors_Tight_ele_1d->Print((plots_dir+"factorized_scalefactors_Tight_ele_1d.png").c_str());
  

  // write to file
  TFile *scalefactors_file = TFile::Open((root_dir+"scalefactors_ele.root").c_str(),"RECREATE");
  unfactorized_scalefactors_Veto_ele   ->Write("unfactorized_scalefactors_Veto_ele"   );
  unfactorized_scalefactors_Loose_ele  ->Write("unfactorized_scalefactors_Loose_ele"  );
  unfactorized_scalefactors_Medium_ele ->Write("unfactorized_scalefactors_Medium_ele" );
  unfactorized_scalefactors_Tight_ele  ->Write("unfactorized_scalefactors_Tight_ele"  );
  factorized_scalefactors_Veto_ele    ->Write("factorized_scalefactors_Veto_ele"  );
  factorized_scalefactors_Loose_ele   ->Write("factorized_scalefactors_Loose_ele" );
  factorized_scalefactors_Medium_ele  ->Write("factorized_scalefactors_Medium_ele");
  factorized_scalefactors_Tight_ele   ->Write("factorized_scalefactors_Tight_ele" );
  factorized_scalefactors_Veto_ele_eta1->Write("factorized_scalefactors_Veto_ele_eta1");
  factorized_scalefactors_Veto_ele_eta2->Write("factorized_scalefactors_Veto_ele_eta2");
  factorized_scalefactors_Loose_ele_eta1->Write("factorized_scalefactors_Loose_ele_eta1");
  factorized_scalefactors_Loose_ele_eta2->Write("factorized_scalefactors_Loose_ele_eta2");
  factorized_scalefactors_Medium_ele_eta1->Write("factorized_scalefactors_Medium_ele_eta1");
  factorized_scalefactors_Medium_ele_eta2->Write("factorized_scalefactors_Medium_ele_eta2");
  factorized_scalefactors_Tight_ele_eta1->Write("factorized_scalefactors_Tight_ele_eta1");
  factorized_scalefactors_Tight_ele_eta2->Write("factorized_scalefactors_Tight_ele_eta2");
  scalefactors_file->Close();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void make_mu_scale_factors(string plots_dir, string root_dir) {
  
  TFile *f_mc   = TFile::Open((root_dir+"DYJetsToLL_efficiencies_muonTnP.root").c_str(),"READ");
  TFile *f_data = TFile::Open((root_dir+"SingleMuon_efficiencies_muonTnP.root").c_str(),"READ");

  // compute Veto scalefactors
  TH2D *eff_Veto_mu_data   = (TH2D*) f_data ->Get("eff_Veto_mu");
  TH2D *eff_Veto_mu_mc     = (TH2D*) f_mc   ->Get("eff_Veto_mu");
  TH2D *unfactorized_scalefactors_Veto_mu = new TH2D(
    "unfactorized_scalefactors_Veto_mu",
    "Muon Veto ID scale factors from Unfactorized efficiencies",
    n_mu_eta_bins, mu_eta_bins,
    n_mu_pt_bins, mu_pt_bins
  );
  (*unfactorized_scalefactors_Veto_mu) = (*eff_Veto_mu_data) / (*eff_Veto_mu_mc);

  // compute Loose scalefactors
  TH2D *eff_Loose_mu_data   = (TH2D*) f_data ->Get("eff_Loose_mu");
  TH2D *eff_Loose_mu_mc     = (TH2D*) f_mc   ->Get("eff_Loose_mu");
  TH2D *unfactorized_scalefactors_Loose_mu = new TH2D(
    "unfactorized_scalefactors_Loose_mu",
    "Muon Loose ID scale factors from Unfactorized efficiencies",
    n_mu_eta_bins, mu_eta_bins,
    n_mu_pt_bins, mu_pt_bins
  );
  (*unfactorized_scalefactors_Loose_mu) = (*eff_Loose_mu_data) / (*eff_Loose_mu_mc);

  // compute Medium scalefactors
  TH2D *eff_Medium_mu_data   = (TH2D*) f_data ->Get("eff_Medium_mu");
  TH2D *eff_Medium_mu_mc     = (TH2D*) f_mc   ->Get("eff_Medium_mu");
  TH2D *unfactorized_scalefactors_Medium_mu = new TH2D(
    "unfactorized_scalefactors_Medium_mu",
    "Muon Medium ID scale factors from Unfactorized efficiencies",
    n_mu_eta_bins, mu_eta_bins,
    n_mu_pt_bins, mu_pt_bins
  );
  (*unfactorized_scalefactors_Medium_mu) = (*eff_Medium_mu_data) / (*eff_Medium_mu_mc);

  // compute Tight scalefactors
  TH2D *eff_Tight_mu_data   = (TH2D*) f_data ->Get("eff_Tight_mu");
  TH2D *eff_Tight_mu_mc     = (TH2D*) f_mc   ->Get("eff_Tight_mu");
  TH2D *unfactorized_scalefactors_Tight_mu = new TH2D(
    "unfactorized_scalefactors_Tight_mu",
    "Muon Tight ID scale factors from Unfactorized efficiencies",
    n_mu_eta_bins, mu_eta_bins,
    n_mu_pt_bins, mu_pt_bins
  );
  (*unfactorized_scalefactors_Tight_mu) = (*eff_Tight_mu_data) / (*eff_Tight_mu_mc);
  
  // Start drawing stuff 
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;
  TLine *oneline=new TLine(10,1,100,1);
  oneline->SetLineColor(1);
  oneline->SetLineStyle(3);
  Int_t mit_red  = 1861; 
  Int_t mit_gray = 1862; 
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);

  // Draw unfactorized scalefactor plots
  TCanvas *c_unfactorized_scalefactors_Veto_mu = new TCanvas("c_unfactorized_scalefactors_Veto_mu","Unfactorized Veto muon scale factors (Data/MC)",800,800);
  unfactorized_scalefactors_Veto_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Veto_mu->SetTitle("Unfactorized Veto muon scale factors (Data/MC)");
  unfactorized_scalefactors_Veto_mu->GetXaxis()->SetTitle("| #eta |");
  unfactorized_scalefactors_Veto_mu->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Veto_mu->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Veto_mu->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Veto_mu->GetYaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Veto_mu->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Veto_mu->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Veto_mu->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Veto_mu->GetYaxis()->SetRangeUser(10,100);
  unfactorized_scalefactors_Veto_mu->SetMinimum(0.7);
  unfactorized_scalefactors_Veto_mu->SetMaximum(1.3);
  unfactorized_scalefactors_Veto_mu->SetMarkerSize(.8);
  //palette_axis = (TPaletteAxis*) unfactorized_scalefactors_Veto_mu->GetListOfFunctions()->FindObject("palette"); 
  //palette_axis->SetLabelSize(0.02);
  c_unfactorized_scalefactors_Veto_mu->Update();
  c_unfactorized_scalefactors_Veto_mu->Print((plots_dir+"unfactorized_scalefactors_Veto_mu.png").c_str());
  TCanvas *c_unfactorized_scalefactors_Loose_mu = new TCanvas("c_unfactorized_scalefactors_Loose_mu","Unfactorized Loose muon scale factors (Data/MC)",800,800);
  unfactorized_scalefactors_Loose_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Loose_mu->SetTitle("Unfactorized Loose muon scale factors (Data/MC)");
  unfactorized_scalefactors_Loose_mu->GetXaxis()->SetTitle("| #eta |");
  unfactorized_scalefactors_Loose_mu->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Loose_mu->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Loose_mu->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Loose_mu->GetYaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Loose_mu->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Loose_mu->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Loose_mu->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Loose_mu->GetYaxis()->SetRangeUser(10,100);
  unfactorized_scalefactors_Loose_mu->SetMinimum(0.7);
  unfactorized_scalefactors_Loose_mu->SetMaximum(1.3);
  unfactorized_scalefactors_Loose_mu->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) unfactorized_scalefactors_Loose_mu->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_unfactorized_scalefactors_Loose_mu->Update();
  c_unfactorized_scalefactors_Loose_mu->Print((plots_dir+"unfactorized_scalefactors_Loose_mu.png").c_str());
  TCanvas *c_unfactorized_scalefactors_Medium_mu = new TCanvas("c_unfactorized_scalefactors_Medium_mu","Unfactorized Medium muon scale factors (Data/MC)",800,800);
  unfactorized_scalefactors_Medium_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Medium_mu->SetTitle("Unfactorized Medium muon scale factors (Data/MC)");
  unfactorized_scalefactors_Medium_mu->GetXaxis()->SetTitle("| #eta |");
  unfactorized_scalefactors_Medium_mu->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Medium_mu->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Medium_mu->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Medium_mu->GetYaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Medium_mu->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Medium_mu->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Medium_mu->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Medium_mu->GetYaxis()->SetRangeUser(10,100);
  unfactorized_scalefactors_Medium_mu->SetMinimum(0.7);
  unfactorized_scalefactors_Medium_mu->SetMaximum(1.3);
  unfactorized_scalefactors_Medium_mu->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) unfactorized_scalefactors_Medium_mu->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_unfactorized_scalefactors_Medium_mu->Update();
  c_unfactorized_scalefactors_Medium_mu->Print((plots_dir+"unfactorized_scalefactors_Medium_mu.png").c_str());
  TCanvas *c_unfactorized_scalefactors_Tight_mu = new TCanvas("c_unfactorized_scalefactors_Tight_mu","Unfactorized Tight muon scale factors (Data/MC)",800,800);
  unfactorized_scalefactors_Tight_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Tight_mu->SetTitle("Unfactorized Tight muon scale factors (Data/MC)");
  unfactorized_scalefactors_Tight_mu->GetXaxis()->SetTitle("| #eta |");
  unfactorized_scalefactors_Tight_mu->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Tight_mu->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Tight_mu->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Tight_mu->GetYaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Tight_mu->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Tight_mu->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Tight_mu->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Tight_mu->GetYaxis()->SetRangeUser(10,100);
  unfactorized_scalefactors_Tight_mu->SetMinimum(0.7);
  unfactorized_scalefactors_Tight_mu->SetMaximum(1.3);
  unfactorized_scalefactors_Tight_mu->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) unfactorized_scalefactors_Tight_mu->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_unfactorized_scalefactors_Tight_mu->Update();
  c_unfactorized_scalefactors_Tight_mu->Print((plots_dir+"unfactorized_scalefactors_Tight_mu.png").c_str());

  TCanvas *c_unfactorized_scalefactors_Veto_mu_1d = new TCanvas("c_unfactorized_scalefactors_Veto_mu_1d","Unfactorized Veto muon scale factors",800,600);
  TH1D *unfactorized_scalefactors_Veto_mu_eta1 = unfactorized_scalefactors_Veto_mu->ProjectionY("unfactorized_scalefactors_Veto_mu_eta1", 1, 1);
  TH1D *unfactorized_scalefactors_Veto_mu_eta2 = unfactorized_scalefactors_Veto_mu->ProjectionY("unfactorized_scalefactors_Veto_mu_eta2", 2, 2);
  //TH1D *unfactorized_scalefactors_Veto_mu_eta3 = unfactorized_scalefactors_Veto_mu->ProjectionY("unfactorized_scalefactors_Veto_mu_eta3", 3, 3);
  unfactorized_scalefactors_Veto_mu_eta1->Draw("P0 E1");
  unfactorized_scalefactors_Veto_mu_eta1->SetTitle("Unfactorized Veto muon scale factors");
  unfactorized_scalefactors_Veto_mu_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Veto_mu_eta1->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Veto_mu_eta1->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Veto_mu_eta1->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Veto_mu_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
  unfactorized_scalefactors_Veto_mu_eta1->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Veto_mu_eta1->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Veto_mu_eta1->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Veto_mu_eta1->SetMinimum(0.8);
  unfactorized_scalefactors_Veto_mu_eta1->SetMaximum(1.2);
  unfactorized_scalefactors_Veto_mu_eta1->SetMarkerStyle(20);
  unfactorized_scalefactors_Veto_mu_eta1->SetMarkerColor(1);
  unfactorized_scalefactors_Veto_mu_eta1->SetLineColor(1);
  unfactorized_scalefactors_Veto_mu_eta2->SetMarkerStyle(21);
  unfactorized_scalefactors_Veto_mu_eta2->SetMarkerColor(1861);
  unfactorized_scalefactors_Veto_mu_eta2->SetLineColor(1861);
  unfactorized_scalefactors_Veto_mu_eta2->Draw("P0 E1 SAME");
  //unfactorized_scalefactors_Veto_mu_eta3->SetMarkerStyle(22);
  //unfactorized_scalefactors_Veto_mu_eta3->SetMarkerColor(1862);
  //unfactorized_scalefactors_Veto_mu_eta3->SetLineColor(1862);
  //unfactorized_scalefactors_Veto_mu_eta3->Draw("P0 E1 SAME");
  oneline->Draw("SAME");
  TLegend *legend_Veto_1d=new TLegend(.6,.7,.85,.85);
  legend_Veto_1d->AddEntry(unfactorized_scalefactors_Veto_mu_eta1,"| #eta | < 1.479", "lp");
  //legend_Veto_1d->AddEntry(unfactorized_scalefactors_Veto_mu_eta2,"1.5 < | #eta | < 2.4", "lp");
  legend_Veto_1d->AddEntry(unfactorized_scalefactors_Veto_mu_eta2,"1.479 < | #eta | < 2.4", "lp");
  //legend_Veto_1d->AddEntry(unfactorized_scalefactors_Veto_mu_eta3,"2.1 < | #eta | < 2.4", "lp");
  legend_Veto_1d->SetFillColor(0);
  legend_Veto_1d->Draw("SAME");
  c_unfactorized_scalefactors_Veto_mu_1d->Update();
  c_unfactorized_scalefactors_Veto_mu_1d->Print((plots_dir+"unfactorized_scalefactors_Veto_mu_1d.png").c_str());
  

  TCanvas *c_unfactorized_scalefactors_Loose_mu_1d = new TCanvas("c_unfactorized_scalefactors_Loose_mu_1d","Unfactorized Loose muon scale factors",800,600);
  TH1D *unfactorized_scalefactors_Loose_mu_eta1 = unfactorized_scalefactors_Loose_mu->ProjectionY("unfactorized_scalefactors_Loose_mu_eta1", 1, 1);
  TH1D *unfactorized_scalefactors_Loose_mu_eta2 = unfactorized_scalefactors_Loose_mu->ProjectionY("unfactorized_scalefactors_Loose_mu_eta2", 2, 2);
  //TH1D *unfactorized_scalefactors_Loose_mu_eta3 = unfactorized_scalefactors_Loose_mu->ProjectionY("unfactorized_scalefactors_Loose_mu_eta3", 3, 3);
  unfactorized_scalefactors_Loose_mu_eta1->Draw("P0 E1");
  unfactorized_scalefactors_Loose_mu_eta1->SetTitle("Unfactorized Loose muon scale factors");
  unfactorized_scalefactors_Loose_mu_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Loose_mu_eta1->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Loose_mu_eta1->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Loose_mu_eta1->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Loose_mu_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
  unfactorized_scalefactors_Loose_mu_eta1->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Loose_mu_eta1->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Loose_mu_eta1->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Loose_mu_eta1->SetMinimum(0.8);
  unfactorized_scalefactors_Loose_mu_eta1->SetMaximum(1.2);
  unfactorized_scalefactors_Loose_mu_eta1->SetMarkerStyle(20);
  unfactorized_scalefactors_Loose_mu_eta1->SetMarkerColor(1);
  unfactorized_scalefactors_Loose_mu_eta1->SetLineColor(1);
  unfactorized_scalefactors_Loose_mu_eta2->SetMarkerStyle(21);
  unfactorized_scalefactors_Loose_mu_eta2->SetMarkerColor(1861);
  unfactorized_scalefactors_Loose_mu_eta2->SetLineColor(1861);
  unfactorized_scalefactors_Loose_mu_eta2->Draw("P0 E1 SAME");
  //unfactorized_scalefactors_Loose_mu_eta3->SetMarkerStyle(22);
  //unfactorized_scalefactors_Loose_mu_eta3->SetMarkerColor(1862);
  //unfactorized_scalefactors_Loose_mu_eta3->SetLineColor(1862);
  //unfactorized_scalefactors_Loose_mu_eta3->Draw("P0 E1 SAME");
  oneline->Draw("SAME");
  TLegend *legend_Loose_1d=new TLegend(.6,.7,.85,.85);
  legend_Loose_1d->AddEntry(unfactorized_scalefactors_Loose_mu_eta1,"| #eta | < 1.479", "lp");
  //legend_Loose_1d->AddEntry(unfactorized_scalefactors_Loose_mu_eta2,"1.5 < | #eta | < 2.4", "lp");
  legend_Loose_1d->AddEntry(unfactorized_scalefactors_Loose_mu_eta2,"1.479 < | #eta | < 2.4", "lp");
  //legend_Loose_1d->AddEntry(unfactorized_scalefactors_Loose_mu_eta3,"2.1 < | #eta | < 2.4", "lp");
  legend_Loose_1d->SetFillColor(0);
  legend_Loose_1d->Draw("SAME");
  c_unfactorized_scalefactors_Loose_mu_1d->Update();
  c_unfactorized_scalefactors_Loose_mu_1d->Print((plots_dir+"unfactorized_scalefactors_Loose_mu_1d.png").c_str());
  

  TCanvas *c_unfactorized_scalefactors_Medium_mu_1d = new TCanvas("c_unfactorized_scalefactors_Medium_mu_1d","Unfactorized Medium muon scale factors",800,600);
  TH1D *unfactorized_scalefactors_Medium_mu_eta1 = unfactorized_scalefactors_Medium_mu->ProjectionY("unfactorized_scalefactors_Medium_mu_eta1", 1, 1);
  TH1D *unfactorized_scalefactors_Medium_mu_eta2 = unfactorized_scalefactors_Medium_mu->ProjectionY("unfactorized_scalefactors_Medium_mu_eta2", 2, 2);
  //TH1D *unfactorized_scalefactors_Medium_mu_eta3 = unfactorized_scalefactors_Medium_mu->ProjectionY("unfactorized_scalefactors_Medium_mu_eta3", 3, 3);
  unfactorized_scalefactors_Medium_mu_eta1->Draw("P0 E1");
  unfactorized_scalefactors_Medium_mu_eta1->SetTitle("Unfactorized Medium muon scale factors");
  unfactorized_scalefactors_Medium_mu_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Medium_mu_eta1->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Medium_mu_eta1->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Medium_mu_eta1->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Medium_mu_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
  unfactorized_scalefactors_Medium_mu_eta1->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Medium_mu_eta1->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Medium_mu_eta1->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Medium_mu_eta1->SetMinimum(0.8);
  unfactorized_scalefactors_Medium_mu_eta1->SetMaximum(1.2);
  unfactorized_scalefactors_Medium_mu_eta1->SetMarkerStyle(20);
  unfactorized_scalefactors_Medium_mu_eta1->SetMarkerColor(1);
  unfactorized_scalefactors_Medium_mu_eta1->SetLineColor(1);
  unfactorized_scalefactors_Medium_mu_eta2->SetMarkerStyle(21);
  unfactorized_scalefactors_Medium_mu_eta2->SetMarkerColor(1861);
  unfactorized_scalefactors_Medium_mu_eta2->SetLineColor(1861);
  unfactorized_scalefactors_Medium_mu_eta2->Draw("P0 E1 SAME");
  //unfactorized_scalefactors_Medium_mu_eta3->SetMarkerStyle(22);
  //unfactorized_scalefactors_Medium_mu_eta3->SetMarkerColor(1862);
  //unfactorized_scalefactors_Medium_mu_eta3->SetLineColor(1862);
  //unfactorized_scalefactors_Medium_mu_eta3->Draw("P0 E1 SAME");
  oneline->Draw("SAME");
  TLegend *legend_Medium_1d=new TLegend(.6,.7,.85,.85);
  legend_Medium_1d->AddEntry(unfactorized_scalefactors_Medium_mu_eta1,"| #eta | < 1.479", "lp");
  //legend_Medium_1d->AddEntry(unfactorized_scalefactors_Medium_mu_eta2,"1.5 < | #eta | < 2.4", "lp");
  legend_Medium_1d->AddEntry(unfactorized_scalefactors_Medium_mu_eta2,"1.479 < | #eta | < 2.4", "lp");
  //legend_Medium_1d->AddEntry(unfactorized_scalefactors_Medium_mu_eta3,"2.1 < | #eta | < 2.4", "lp");
  legend_Medium_1d->SetFillColor(0);
  legend_Medium_1d->Draw("SAME");
  c_unfactorized_scalefactors_Medium_mu_1d->Update();
  c_unfactorized_scalefactors_Medium_mu_1d->Print((plots_dir+"unfactorized_scalefactors_Medium_mu_1d.png").c_str());
  

  TCanvas *c_unfactorized_scalefactors_Tight_mu_1d = new TCanvas("c_unfactorized_scalefactors_Tight_mu_1d","Unfactorized Tight muon scale factors",800,600);
  TH1D *unfactorized_scalefactors_Tight_mu_eta1 = unfactorized_scalefactors_Tight_mu->ProjectionY("unfactorized_scalefactors_Tight_mu_eta1", 1, 1);
  TH1D *unfactorized_scalefactors_Tight_mu_eta2 = unfactorized_scalefactors_Tight_mu->ProjectionY("unfactorized_scalefactors_Tight_mu_eta2", 2, 2);
  //TH1D *unfactorized_scalefactors_Tight_mu_eta3 = unfactorized_scalefactors_Tight_mu->ProjectionY("unfactorized_scalefactors_Tight_mu_eta3", 3, 3);
  unfactorized_scalefactors_Tight_mu_eta1->Draw("P0 E1");
  unfactorized_scalefactors_Tight_mu_eta1->SetTitle("Unfactorized Tight muon scale factors");
  unfactorized_scalefactors_Tight_mu_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
  unfactorized_scalefactors_Tight_mu_eta1->GetXaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Tight_mu_eta1->GetXaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Tight_mu_eta1->GetXaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Tight_mu_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
  unfactorized_scalefactors_Tight_mu_eta1->GetYaxis()->SetTitleOffset(0.9);
  unfactorized_scalefactors_Tight_mu_eta1->GetYaxis()->SetTitleSize(0.04);
  unfactorized_scalefactors_Tight_mu_eta1->GetYaxis()->SetLabelSize(0.02);
  unfactorized_scalefactors_Tight_mu_eta1->SetMinimum(0.8);
  unfactorized_scalefactors_Tight_mu_eta1->SetMaximum(1.2);
  unfactorized_scalefactors_Tight_mu_eta1->SetMarkerStyle(20);
  unfactorized_scalefactors_Tight_mu_eta1->SetMarkerColor(1);
  unfactorized_scalefactors_Tight_mu_eta1->SetLineColor(1);
  unfactorized_scalefactors_Tight_mu_eta2->SetMarkerStyle(21);
  unfactorized_scalefactors_Tight_mu_eta2->SetMarkerColor(1861);
  unfactorized_scalefactors_Tight_mu_eta2->SetLineColor(1861);
  unfactorized_scalefactors_Tight_mu_eta2->Draw("P0 E1 SAME");
  //unfactorized_scalefactors_Tight_mu_eta3->SetMarkerStyle(22);
  //unfactorized_scalefactors_Tight_mu_eta3->SetMarkerColor(1862);
  //unfactorized_scalefactors_Tight_mu_eta3->SetLineColor(1862);
  //unfactorized_scalefactors_Tight_mu_eta3->Draw("P0 E1 SAME");
  oneline->Draw("SAME");
  TLegend *legend_Tight_1d=new TLegend(.6,.7,.85,.85);
  legend_Tight_1d->AddEntry(unfactorized_scalefactors_Tight_mu_eta1,"| #eta | < 1.479", "lp");
  //legend_Tight_1d->AddEntry(unfactorized_scalefactors_Tight_mu_eta2,"1.5 < | #eta | < 2.4", "lp");
  legend_Tight_1d->AddEntry(unfactorized_scalefactors_Tight_mu_eta2,"1.479 < | #eta | < 2.4", "lp");
  //legend_Tight_1d->AddEntry(unfactorized_scalefactors_Tight_mu_eta3,"2.1 < | #eta | < 2.4", "lp");
  legend_Tight_1d->SetFillColor(0);
  legend_Tight_1d->Draw("SAME");
  c_unfactorized_scalefactors_Tight_mu_1d->Update();
  c_unfactorized_scalefactors_Tight_mu_1d->Print((plots_dir+"unfactorized_scalefactors_Tight_mu_1d.png").c_str());
  
  
  TFile *scalefactors_file = TFile::Open((root_dir+"scalefactors_mu.root").c_str(),"RECREATE");
  unfactorized_scalefactors_Veto_mu   ->Write("unfactorized_scalefactors_Veto_mu"  );
  unfactorized_scalefactors_Loose_mu  ->Write("unfactorized_scalefactors_Loose_mu" );
  unfactorized_scalefactors_Medium_mu ->Write("unfactorized_scalefactors_Medium_mu");
  unfactorized_scalefactors_Tight_mu  ->Write("unfactorized_scalefactors_Tight_mu" );
  unfactorized_scalefactors_Veto_mu_eta1->Write("unfactorized_scalefactors_Veto_mu_eta1");
  unfactorized_scalefactors_Veto_mu_eta2->Write("unfactorized_scalefactors_Veto_mu_eta2");
//  unfactorized_scalefactors_Veto_mu_eta3->Write("unfactorized_scalefactors_Veto_mu_eta3");
  unfactorized_scalefactors_Loose_mu_eta1->Write("unfactorized_scalefactors_Loose_mu_eta1");
  unfactorized_scalefactors_Loose_mu_eta2->Write("unfactorized_scalefactors_Loose_mu_eta2");
//  unfactorized_scalefactors_Loose_mu_eta3->Write("unfactorized_scalefactors_Loose_mu_eta3");
  unfactorized_scalefactors_Medium_mu_eta1->Write("unfactorized_scalefactors_Medium_mu_eta1");
  unfactorized_scalefactors_Medium_mu_eta2->Write("unfactorized_scalefactors_Medium_mu_eta2");
//  unfactorized_scalefactors_Medium_mu_eta3->Write("unfactorized_scalefactors_Medium_mu_eta3");
  unfactorized_scalefactors_Tight_mu_eta1->Write("unfactorized_scalefactors_Tight_mu_eta1");
  unfactorized_scalefactors_Tight_mu_eta2->Write("unfactorized_scalefactors_Tight_mu_eta2");
//  unfactorized_scalefactors_Tight_mu_eta3->Write("unfactorized_scalefactors_Tight_mu_eta3");
  scalefactors_file->Close();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// estimate systematic uncertainty by taking absolute difference in scale factors from two methods
void estimate_systematics(string plots_dir_method1, string plots_dir_method2, string root_dir_method1, string root_dir_method2) {
  
  TFile *ele_scalefactors_file_method1 = TFile::Open((root_dir_method1+"scalefactors_ele.root").c_str(),"UPDATE");
  TFile *ele_scalefactors_file_method2 = TFile::Open((root_dir_method2+"scalefactors_ele.root").c_str(),"UPDATE");
  TFile *mu_scalefactors_file_method1 = TFile::Open((root_dir_method1+"scalefactors_mu.root").c_str(),"UPDATE");
  TFile *mu_scalefactors_file_method2 = TFile::Open((root_dir_method2+"scalefactors_mu.root").c_str(),"UPDATE");
  if(!ele_scalefactors_file_method1 || !ele_scalefactors_file_method2 || !mu_scalefactors_file_method1 || !mu_scalefactors_file_method2) {
    printf("error opening scalefactors files\n"); return;
  }

  
  TH2D *unfactorized_scalefactors_Veto_mu_method1   = (TH2D*)mu_scalefactors_file_method1->Get("unfactorized_scalefactors_Veto_mu");
  TH2D *unfactorized_scalefactors_Loose_mu_method1  = (TH2D*)mu_scalefactors_file_method1->Get("unfactorized_scalefactors_Loose_mu");
  TH2D *unfactorized_scalefactors_Medium_mu_method1 = (TH2D*)mu_scalefactors_file_method1->Get("unfactorized_scalefactors_Medium_mu");
  TH2D *unfactorized_scalefactors_Tight_mu_method1  = (TH2D*)mu_scalefactors_file_method1->Get("unfactorized_scalefactors_Tight_mu");
  TH2D *unfactorized_scalefactors_Veto_mu_method2   = (TH2D*)mu_scalefactors_file_method2->Get("unfactorized_scalefactors_Veto_mu");
  TH2D *unfactorized_scalefactors_Loose_mu_method2  = (TH2D*)mu_scalefactors_file_method2->Get("unfactorized_scalefactors_Loose_mu");
  TH2D *unfactorized_scalefactors_Medium_mu_method2 = (TH2D*)mu_scalefactors_file_method2->Get("unfactorized_scalefactors_Medium_mu");
  TH2D *unfactorized_scalefactors_Tight_mu_method2  = (TH2D*)mu_scalefactors_file_method2->Get("unfactorized_scalefactors_Tight_mu");
  TH2D *factorized_scalefactors_Veto_ele_method1   = (TH2D*)ele_scalefactors_file_method1->Get("factorized_scalefactors_Veto_ele");
  TH2D *factorized_scalefactors_Loose_ele_method1  = (TH2D*)ele_scalefactors_file_method1->Get("factorized_scalefactors_Loose_ele");
  TH2D *factorized_scalefactors_Medium_ele_method1 = (TH2D*)ele_scalefactors_file_method1->Get("factorized_scalefactors_Medium_ele");
  TH2D *factorized_scalefactors_Tight_ele_method1  = (TH2D*)ele_scalefactors_file_method1->Get("factorized_scalefactors_Tight_ele");
  TH2D *factorized_scalefactors_Veto_ele_method2   = (TH2D*)ele_scalefactors_file_method2->Get("factorized_scalefactors_Veto_ele");
  TH2D *factorized_scalefactors_Loose_ele_method2  = (TH2D*)ele_scalefactors_file_method2->Get("factorized_scalefactors_Loose_ele");
  TH2D *factorized_scalefactors_Medium_ele_method2 = (TH2D*)ele_scalefactors_file_method2->Get("factorized_scalefactors_Medium_ele");
  TH2D *factorized_scalefactors_Tight_ele_method2  = (TH2D*)ele_scalefactors_file_method2->Get("factorized_scalefactors_Tight_ele");
  
  TH2D *absDiff_unfactorized_scalefactors_Veto_mu = new TH2D(
    "absDiff_unfactorized_scalefactors_Veto_mu",
    "Absolute difference between unfactorized Muon Veto scalefactors computed 2 ways",
    n_mu_eta_bins, mu_eta_bins,
    n_mu_pt_bins, mu_pt_bins
  );
  TH2D *absDiff_unfactorized_scalefactors_Loose_mu = new TH2D(
    "absDiff_unfactorized_scalefactors_Loose_mu",
    "Absolute difference between unfactorized Muon Loose scalefactors computed 2 ways",
    n_mu_eta_bins, mu_eta_bins,
    n_mu_pt_bins, mu_pt_bins
  );
  TH2D *absDiff_unfactorized_scalefactors_Medium_mu = new TH2D(
    "absDiff_unfactorized_scalefactors_Medium_mu",
    "Absolute difference between unfactorized Muon Medium scalefactors computed 2 ways",
    n_mu_eta_bins, mu_eta_bins,
    n_mu_pt_bins, mu_pt_bins
  );
  TH2D *absDiff_unfactorized_scalefactors_Tight_mu = new TH2D(
    "absDiff_unfactorized_scalefactors_Tight_mu",
    "Absolute difference between unfactorized Muon Tight scalefactors computed 2 ways",
    n_mu_eta_bins, mu_eta_bins,
    n_mu_pt_bins, mu_pt_bins
  );
  TH2D *absDiff_factorized_scalefactors_Veto_ele = new TH2D(
    "absDiff_factorized_scalefactors_Veto_ele",
    "Absolute difference between factorized Electron Veto scalefactors computed 2 ways",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  TH2D *absDiff_factorized_scalefactors_Loose_ele = new TH2D(
    "absDiff_factorized_scalefactors_Loose_ele",
    "Absolute difference between factorized Electron Loose scalefactors computed 2 ways",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  TH2D *absDiff_factorized_scalefactors_Medium_ele = new TH2D(
    "absDiff_factorized_scalefactors_Medium_ele",
    "Absolute difference between factorized Electron Medium scalefactors computed 2 ways",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  TH2D *absDiff_factorized_scalefactors_Tight_ele = new TH2D(
    "absDiff_factorized_scalefactors_Tight_ele",
    "Absolute difference between factorized Electron Tight scalefactors computed 2 ways",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );

  // compute difference
  (*absDiff_unfactorized_scalefactors_Veto_mu)   = (*unfactorized_scalefactors_Veto_mu_method1)   - (*unfactorized_scalefactors_Veto_mu_method2);
  (*absDiff_unfactorized_scalefactors_Loose_mu)  = (*unfactorized_scalefactors_Loose_mu_method1)  - (*unfactorized_scalefactors_Loose_mu_method2);
  (*absDiff_unfactorized_scalefactors_Medium_mu) = (*unfactorized_scalefactors_Medium_mu_method1) - (*unfactorized_scalefactors_Medium_mu_method2);
  (*absDiff_unfactorized_scalefactors_Tight_mu)  = (*unfactorized_scalefactors_Tight_mu_method1)  - (*unfactorized_scalefactors_Tight_mu_method2);
  (*absDiff_factorized_scalefactors_Veto_ele)   = (*factorized_scalefactors_Veto_ele_method1)   - (*factorized_scalefactors_Veto_ele_method2);
  (*absDiff_factorized_scalefactors_Loose_ele)  = (*factorized_scalefactors_Loose_ele_method1)  - (*factorized_scalefactors_Loose_ele_method2);
  (*absDiff_factorized_scalefactors_Medium_ele) = (*factorized_scalefactors_Medium_ele_method1) - (*factorized_scalefactors_Medium_ele_method2);
  (*absDiff_factorized_scalefactors_Tight_ele)  = (*factorized_scalefactors_Tight_ele_method1)  - (*factorized_scalefactors_Tight_ele_method2);
  // take its absolute value
  for(int i = 1; i<=n_mu_eta_bins; i++) { for(int j = 1; j<=n_mu_pt_bins; j++) {
    int nbin = absDiff_unfactorized_scalefactors_Veto_mu->GetBin(i,j);
    absDiff_unfactorized_scalefactors_Veto_mu   ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Veto_mu  ->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Loose_mu  ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Loose_mu ->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Medium_mu ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Medium_mu->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Tight_mu  ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Tight_mu ->GetBinContent(nbin)));  
  }}
  for(int i = 1; i<=n_ele_eta_bins; i++) { for(int j = 1; j<=n_ele_pt_bins; j++) {
    int nbin = absDiff_factorized_scalefactors_Veto_ele->GetBin(i,j);
    absDiff_factorized_scalefactors_Veto_ele   ->SetBinContent(nbin, TMath::Abs(absDiff_factorized_scalefactors_Veto_ele  ->GetBinContent(nbin)));  
    absDiff_factorized_scalefactors_Loose_ele  ->SetBinContent(nbin, TMath::Abs(absDiff_factorized_scalefactors_Loose_ele ->GetBinContent(nbin)));  
    absDiff_factorized_scalefactors_Medium_ele ->SetBinContent(nbin, TMath::Abs(absDiff_factorized_scalefactors_Medium_ele->GetBinContent(nbin)));  
    absDiff_factorized_scalefactors_Tight_ele  ->SetBinContent(nbin, TMath::Abs(absDiff_factorized_scalefactors_Tight_ele ->GetBinContent(nbin)));  
  }}

  // Start drawing stuff 
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;
  TLine *oneline=new TLine(10,1,100,1);
  oneline->SetLineColor(1);
  oneline->SetLineStyle(3);
  Int_t mit_red  = 1861; 
  Int_t mit_gray = 1862; 
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);

  // Draw absolute difference plots
  TCanvas *c_absDiff_unfactorized_scalefactors_Veto_mu = new TCanvas("c_absDiff_unfactorized_scalefactors_Veto_mu","Method |difference| for Veto muon scale factors",800,800);
  absDiff_unfactorized_scalefactors_Veto_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_unfactorized_scalefactors_Veto_mu->SetTitle("Method |difference| for Veto Muon scale factors");
  absDiff_unfactorized_scalefactors_Veto_mu->GetXaxis()->SetTitle("| #eta |");
  absDiff_unfactorized_scalefactors_Veto_mu->GetXaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Veto_mu->GetXaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Veto_mu->GetXaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Veto_mu->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_unfactorized_scalefactors_Veto_mu->GetYaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Veto_mu->GetYaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Veto_mu->GetYaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Veto_mu->GetYaxis()->SetRangeUser(10,100);
  absDiff_unfactorized_scalefactors_Veto_mu->SetMinimum(0);
  absDiff_unfactorized_scalefactors_Veto_mu->SetMaximum(0.2);
  absDiff_unfactorized_scalefactors_Veto_mu->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_unfactorized_scalefactors_Veto_mu->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_unfactorized_scalefactors_Veto_mu->Update();
  c_absDiff_unfactorized_scalefactors_Veto_mu->Print((plots_dir_method1+"absDiff_unfactorized_scalefactors_Veto_mu.png").c_str());
  c_absDiff_unfactorized_scalefactors_Veto_mu->Print((plots_dir_method2+"absDiff_unfactorized_scalefactors_Veto_mu.png").c_str());

  TCanvas *c_absDiff_unfactorized_scalefactors_Loose_mu = new TCanvas("c_absDiff_unfactorized_scalefactors_Loose_mu","Method |difference| for Loose muon scale factors",800,800);
  absDiff_unfactorized_scalefactors_Loose_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_unfactorized_scalefactors_Loose_mu->SetTitle("Method |difference| for Loose Muon scale factors");
  absDiff_unfactorized_scalefactors_Loose_mu->GetXaxis()->SetTitle("| #eta |");
  absDiff_unfactorized_scalefactors_Loose_mu->GetXaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Loose_mu->GetXaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Loose_mu->GetXaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Loose_mu->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_unfactorized_scalefactors_Loose_mu->GetYaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Loose_mu->GetYaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Loose_mu->GetYaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Loose_mu->GetYaxis()->SetRangeUser(10,100);
  absDiff_unfactorized_scalefactors_Loose_mu->SetMinimum(0);
  absDiff_unfactorized_scalefactors_Loose_mu->SetMaximum(0.2);
  absDiff_unfactorized_scalefactors_Loose_mu->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_unfactorized_scalefactors_Loose_mu->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_unfactorized_scalefactors_Loose_mu->Update();
  c_absDiff_unfactorized_scalefactors_Loose_mu->Print((plots_dir_method1+"absDiff_unfactorized_scalefactors_Loose_mu.png").c_str());
  c_absDiff_unfactorized_scalefactors_Loose_mu->Print((plots_dir_method2+"absDiff_unfactorized_scalefactors_Loose_mu.png").c_str());

  TCanvas *c_absDiff_unfactorized_scalefactors_Medium_mu = new TCanvas("c_absDiff_unfactorized_scalefactors_Medium_mu","Method |difference| for Medium muon scale factors",800,800);
  absDiff_unfactorized_scalefactors_Medium_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_unfactorized_scalefactors_Medium_mu->SetTitle("Method |difference| for Medium Muon scale factors");
  absDiff_unfactorized_scalefactors_Medium_mu->GetXaxis()->SetTitle("| #eta |");
  absDiff_unfactorized_scalefactors_Medium_mu->GetXaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Medium_mu->GetXaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Medium_mu->GetXaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Medium_mu->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_unfactorized_scalefactors_Medium_mu->GetYaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Medium_mu->GetYaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Medium_mu->GetYaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Medium_mu->GetYaxis()->SetRangeUser(10,100);
  absDiff_unfactorized_scalefactors_Medium_mu->SetMinimum(0);
  absDiff_unfactorized_scalefactors_Medium_mu->SetMaximum(0.2);
  absDiff_unfactorized_scalefactors_Medium_mu->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_unfactorized_scalefactors_Medium_mu->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_unfactorized_scalefactors_Medium_mu->Update();
  c_absDiff_unfactorized_scalefactors_Medium_mu->Print((plots_dir_method1+"absDiff_unfactorized_scalefactors_Medium_mu.png").c_str());
  c_absDiff_unfactorized_scalefactors_Medium_mu->Print((plots_dir_method2+"absDiff_unfactorized_scalefactors_Medium_mu.png").c_str());

  TCanvas *c_absDiff_unfactorized_scalefactors_Tight_mu = new TCanvas("c_absDiff_unfactorized_scalefactors_Tight_mu","Method |difference| for Tight muon scale factors",800,800);
  absDiff_unfactorized_scalefactors_Tight_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_unfactorized_scalefactors_Tight_mu->SetTitle("Method |difference| for Tight Muon scale factors");
  absDiff_unfactorized_scalefactors_Tight_mu->GetXaxis()->SetTitle("| #eta |");
  absDiff_unfactorized_scalefactors_Tight_mu->GetXaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Tight_mu->GetXaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Tight_mu->GetXaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Tight_mu->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_unfactorized_scalefactors_Tight_mu->GetYaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Tight_mu->GetYaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Tight_mu->GetYaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Tight_mu->GetYaxis()->SetRangeUser(10,100);
  absDiff_unfactorized_scalefactors_Tight_mu->SetMinimum(0);
  absDiff_unfactorized_scalefactors_Tight_mu->SetMaximum(0.2);
  absDiff_unfactorized_scalefactors_Tight_mu->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_unfactorized_scalefactors_Tight_mu->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_unfactorized_scalefactors_Tight_mu->Update();
  c_absDiff_unfactorized_scalefactors_Tight_mu->Print((plots_dir_method1+"absDiff_unfactorized_scalefactors_Tight_mu.png").c_str());
  c_absDiff_unfactorized_scalefactors_Tight_mu->Print((plots_dir_method2+"absDiff_unfactorized_scalefactors_Tight_mu.png").c_str());

  TCanvas *c_absDiff_factorized_scalefactors_Veto_ele = new TCanvas("c_absDiff_factorized_scalefactors_Veto_ele","Method |difference| for Veto electron scale factors",800,800);
  absDiff_factorized_scalefactors_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_factorized_scalefactors_Veto_ele->SetTitle("Method |difference| for Veto Electron scale factors");
  absDiff_factorized_scalefactors_Veto_ele->GetXaxis()->SetTitle("| #eta |");
  absDiff_factorized_scalefactors_Veto_ele->GetXaxis()->SetTitleOffset(0.9);
  absDiff_factorized_scalefactors_Veto_ele->GetXaxis()->SetTitleSize(0.04);
  absDiff_factorized_scalefactors_Veto_ele->GetXaxis()->SetLabelSize(0.02);
  absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetTitleOffset(0.9);
  absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetTitleSize(0.04);
  absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetLabelSize(0.02);
  absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetRangeUser(10,100);
  absDiff_factorized_scalefactors_Veto_ele->SetMinimum(0);
  absDiff_factorized_scalefactors_Veto_ele->SetMaximum(0.2);
  absDiff_factorized_scalefactors_Veto_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_factorized_scalefactors_Veto_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_factorized_scalefactors_Veto_ele->Update();
  c_absDiff_factorized_scalefactors_Veto_ele->Print((plots_dir_method1+"absDiff_factorized_scalefactors_Veto_ele.png").c_str());
  c_absDiff_factorized_scalefactors_Veto_ele->Print((plots_dir_method2+"absDiff_factorized_scalefactors_Veto_ele.png").c_str());

  TCanvas *c_absDiff_factorized_scalefactors_Loose_ele = new TCanvas("c_absDiff_factorized_scalefactors_Loose_ele","Method |difference| for Loose electron scale factors",800,800);
  absDiff_factorized_scalefactors_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_factorized_scalefactors_Loose_ele->SetTitle("Method |difference| for Loose Electron scale factors");
  absDiff_factorized_scalefactors_Loose_ele->GetXaxis()->SetTitle("| #eta |");
  absDiff_factorized_scalefactors_Loose_ele->GetXaxis()->SetTitleOffset(0.9);
  absDiff_factorized_scalefactors_Loose_ele->GetXaxis()->SetTitleSize(0.04);
  absDiff_factorized_scalefactors_Loose_ele->GetXaxis()->SetLabelSize(0.02);
  absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetTitleOffset(0.9);
  absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetTitleSize(0.04);
  absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetLabelSize(0.02);
  absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetRangeUser(10,100);
  absDiff_factorized_scalefactors_Loose_ele->SetMinimum(0);
  absDiff_factorized_scalefactors_Loose_ele->SetMaximum(0.2);
  absDiff_factorized_scalefactors_Loose_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_factorized_scalefactors_Loose_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_factorized_scalefactors_Loose_ele->Update();
  c_absDiff_factorized_scalefactors_Loose_ele->Print((plots_dir_method1+"absDiff_factorized_scalefactors_Loose_ele.png").c_str());
  c_absDiff_factorized_scalefactors_Loose_ele->Print((plots_dir_method2+"absDiff_factorized_scalefactors_Loose_ele.png").c_str());

  TCanvas *c_absDiff_factorized_scalefactors_Medium_ele = new TCanvas("c_absDiff_factorized_scalefactors_Medium_ele","Method |difference| for Medium electron scale factors",800,800);
  absDiff_factorized_scalefactors_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_factorized_scalefactors_Medium_ele->SetTitle("Method |difference| for Medium Electron scale factors");
  absDiff_factorized_scalefactors_Medium_ele->GetXaxis()->SetTitle("| #eta |");
  absDiff_factorized_scalefactors_Medium_ele->GetXaxis()->SetTitleOffset(0.9);
  absDiff_factorized_scalefactors_Medium_ele->GetXaxis()->SetTitleSize(0.04);
  absDiff_factorized_scalefactors_Medium_ele->GetXaxis()->SetLabelSize(0.02);
  absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetTitleOffset(0.9);
  absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetTitleSize(0.04);
  absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetLabelSize(0.02);
  absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetRangeUser(10,100);
  absDiff_factorized_scalefactors_Medium_ele->SetMinimum(0);
  absDiff_factorized_scalefactors_Medium_ele->SetMaximum(0.2);
  absDiff_factorized_scalefactors_Medium_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_factorized_scalefactors_Medium_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_factorized_scalefactors_Medium_ele->Update();
  c_absDiff_factorized_scalefactors_Medium_ele->Print((plots_dir_method1+"absDiff_factorized_scalefactors_Medium_ele.png").c_str());
  c_absDiff_factorized_scalefactors_Medium_ele->Print((plots_dir_method2+"absDiff_factorized_scalefactors_Medium_ele.png").c_str());

  TCanvas *c_absDiff_factorized_scalefactors_Tight_ele = new TCanvas("c_absDiff_factorized_scalefactors_Tight_ele","Method |difference| for Tight electron scale factors",800,800);
  absDiff_factorized_scalefactors_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_factorized_scalefactors_Tight_ele->SetTitle("Method |difference| for Tight Electron scale factors");
  absDiff_factorized_scalefactors_Tight_ele->GetXaxis()->SetTitle("| #eta |");
  absDiff_factorized_scalefactors_Tight_ele->GetXaxis()->SetTitleOffset(0.9);
  absDiff_factorized_scalefactors_Tight_ele->GetXaxis()->SetTitleSize(0.04);
  absDiff_factorized_scalefactors_Tight_ele->GetXaxis()->SetLabelSize(0.02);
  absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetTitleOffset(0.9);
  absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetTitleSize(0.04);
  absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetLabelSize(0.02);
  absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetRangeUser(10,100);
  absDiff_factorized_scalefactors_Tight_ele->SetMinimum(0);
  absDiff_factorized_scalefactors_Tight_ele->SetMaximum(0.2);
  absDiff_factorized_scalefactors_Tight_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_factorized_scalefactors_Tight_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_factorized_scalefactors_Tight_ele->Update();
  c_absDiff_factorized_scalefactors_Tight_ele->Print((plots_dir_method1+"absDiff_factorized_scalefactors_Tight_ele.png").c_str());
  c_absDiff_factorized_scalefactors_Tight_ele->Print((plots_dir_method2+"absDiff_factorized_scalefactors_Tight_ele.png").c_str());

  mu_scalefactors_file_method1->cd(); 
  mu_scalefactors_file_method1->Delete("absDiff_unfactorized_scalefactors_Veto_mu;*");
  mu_scalefactors_file_method1->Delete("absDiff_unfactorized_scalefactors_Loose_mu;*");
  mu_scalefactors_file_method1->Delete("absDiff_unfactorized_scalefactors_Medium_mu;*");
  mu_scalefactors_file_method1->Delete("absDiff_unfactorized_scalefactors_Tight_mu;*");
  absDiff_unfactorized_scalefactors_Veto_mu->Write("absDiff_unfactorized_scalefactors_Veto_mu");
  absDiff_unfactorized_scalefactors_Loose_mu->Write("absDiff_unfactorized_scalefactors_Loose_mu");
  absDiff_unfactorized_scalefactors_Medium_mu->Write("absDiff_unfactorized_scalefactors_Medium_mu");
  absDiff_unfactorized_scalefactors_Tight_mu->Write("absDiff_unfactorized_scalefactors_Tight_mu");
  mu_scalefactors_file_method2->cd();
  mu_scalefactors_file_method2->Delete("absDiff_unfactorized_scalefactors_Veto_mu;*");
  mu_scalefactors_file_method2->Delete("absDiff_unfactorized_scalefactors_Loose_mu;*");
  mu_scalefactors_file_method2->Delete("absDiff_unfactorized_scalefactors_Medium_mu;*");
  mu_scalefactors_file_method2->Delete("absDiff_unfactorized_scalefactors_Tight_mu;*");
  absDiff_unfactorized_scalefactors_Veto_mu->Write("absDiff_unfactorized_scalefactors_Veto_mu");
  absDiff_unfactorized_scalefactors_Loose_mu->Write("absDiff_unfactorized_scalefactors_Loose_mu");
  absDiff_unfactorized_scalefactors_Medium_mu->Write("absDiff_unfactorized_scalefactors_Medium_mu");
  absDiff_unfactorized_scalefactors_Tight_mu->Write("absDiff_unfactorized_scalefactors_Tight_mu");
  ele_scalefactors_file_method1->cd();
  ele_scalefactors_file_method1->Delete("absDiff_factorized_scalefactors_Veto_ele;*");
  ele_scalefactors_file_method1->Delete("absDiff_factorized_scalefactors_Loose_ele;*");
  ele_scalefactors_file_method1->Delete("absDiff_factorized_scalefactors_Medium_ele;*");
  ele_scalefactors_file_method1->Delete("absDiff_factorized_scalefactors_Tight_ele;*");
  absDiff_factorized_scalefactors_Veto_ele->Write("absDiff_factorized_scalefactors_Veto_ele");
  absDiff_factorized_scalefactors_Loose_ele->Write("absDiff_factorized_scalefactors_Loose_ele");
  absDiff_factorized_scalefactors_Medium_ele->Write("absDiff_factorized_scalefactors_Medium_ele");
  absDiff_factorized_scalefactors_Tight_ele->Write("absDiff_factorized_scalefactors_Tight_ele");
  ele_scalefactors_file_method2->cd();
  ele_scalefactors_file_method2->Delete("absDiff_factorized_scalefactors_Veto_ele;*");
  ele_scalefactors_file_method2->Delete("absDiff_factorized_scalefactors_Loose_ele;*");
  ele_scalefactors_file_method2->Delete("absDiff_factorized_scalefactors_Medium_ele;*");
  ele_scalefactors_file_method2->Delete("absDiff_factorized_scalefactors_Tight_ele;*");
  absDiff_factorized_scalefactors_Veto_ele->Write("absDiff_factorized_scalefactors_Veto_ele");
  absDiff_factorized_scalefactors_Loose_ele->Write("absDiff_factorized_scalefactors_Loose_ele");
  absDiff_factorized_scalefactors_Medium_ele->Write("absDiff_factorized_scalefactors_Medium_ele");
  absDiff_factorized_scalefactors_Tight_ele->Write("absDiff_factorized_scalefactors_Tight_ele");
  printf("* press any key to continue * \n\n");
  std::cin.ignore();

  ele_scalefactors_file_method1->Close();
  ele_scalefactors_file_method2->Close();
  mu_scalefactors_file_method1->Close();
  mu_scalefactors_file_method2->Close();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void wholething(string plots_dir, string root_dir) {
  efficiencies_mc_mu(plots_dir, root_dir);
  efficiencies_data_mu(plots_dir, root_dir);
  factorize_mc_ele(plots_dir, root_dir);
  factorize_data_ele(plots_dir, root_dir);
  make_ele_scale_factors(plots_dir, root_dir);
  make_mu_scale_factors(plots_dir, root_dir);
}
