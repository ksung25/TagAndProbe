#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TPaletteAxis.h>

void mitPalette()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;
   Double_t Red[3]    = { 1, 138./255., 163/255.};
   Double_t Green[3]  = { 1, 139./255., 31/255.};
   Double_t Blue[3]   = { 1, 140./255., 52/255.};
   Double_t Length[3] = { 0.00, 0.35, 1.00 };
   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetPalette(50,colors);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void factorize_mc_ele(string plots_dir, string root_dir) {
  Float_t ele_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t ele_eta_bins[] = {0, 1.479, 2.4};
  Int_t n_ele_pt_bins=12;
  Int_t n_ele_eta_bins=2;
  Float_t mu_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t mu_eta_bins[] = {0, 1.5, 2.1, 2.4};
  Int_t n_mu_pt_bins=12;
  Int_t n_mu_eta_bins=3;
  
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
  TCanvas *c_absdiff_Veto_ele = new TCanvas("c_absdiff_Veto_ele","Difference from factorized efficiency, Electron Veto ID (MC)",800,800);
  absdiff_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Veto_ele->SetTitle("Difference from factorized efficiency, Electron Veto ID (MC)");
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
  TCanvas *c_absdiff_Loose_ele = new TCanvas("c_absdiff_Loose_ele","Difference from factorized efficiency, Electron Loose ID (MC)",800,800);
  absdiff_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Loose_ele->SetTitle("Difference from factorized efficiency, Electron Loose ID (MC)");
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
  TCanvas *c_absdiff_Medium_ele = new TCanvas("c_absdiff_Medium_ele","Difference from factorized efficiency, Electron Medium ID (MC)",800,800);
  absdiff_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Medium_ele->SetTitle("Difference from factorized efficiency, Electron Medium ID (MC)");
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
  TCanvas *c_absdiff_Tight_ele = new TCanvas("c_absdiff_Tight_ele","Difference from factorized efficiency, Electron Tight ID (MC)",800,800);
  absdiff_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Tight_ele->SetTitle("Difference from factorized efficiency, Electron Tight ID (MC)");
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
  TCanvas *c_eff_Veto_ele = new TCanvas("c_eff_Veto_ele","Unfactorized efficiency, Electron Veto ID (MC)",800,800);
  eff_BaselineToVeto_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToVeto_electronTnP->SetTitle("Unfactorized efficiency, Electron Veto ID (MC)");
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
  TCanvas *c_eff_Loose_ele = new TCanvas("c_eff_Loose_ele","Unfactorized efficiency, Electron Loose ID (MC)",800,800);
  eff_BaselineToLoose_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToLoose_electronTnP->SetTitle("Unfactorized efficiency, Electron Loose ID (MC)");
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
  TCanvas *c_eff_Medium_ele = new TCanvas("c_eff_Medium_ele","Unfactorized efficiency, Electron Medium ID (MC)",800,800);
  eff_BaselineToMedium_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToMedium_electronTnP->SetTitle("Unfactorized efficiency, Electron Medium ID (MC)");
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
  TCanvas *c_eff_Tight_ele = new TCanvas("c_eff_Tight_ele","Unfactorized efficiency, Electron Tight ID (MC)",800,800);
  eff_BaselineToTight_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToTight_electronTnP->SetTitle("Unfactorized efficiency, Electron Tight ID (MC)");
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
  TCanvas *c_factEff_Veto_ele = new TCanvas("c_factEff_Veto_ele","Factorized efficiency, Electron Veto ID (MC)",800,800);
  factEff_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Veto_ele->SetTitle("Factorized efficiency, Electron Veto ID (MC)");
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
  TCanvas *c_factEff_Loose_ele = new TCanvas("c_factEff_Loose_ele","Factorized efficiency, Electron Loose ID (MC)",800,800);
  factEff_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Loose_ele->SetTitle("Factorized efficiency, Electron Loose ID (MC)");
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
  TCanvas *c_factEff_Medium_ele = new TCanvas("c_factEff_Medium_ele","Factorized efficiency, Electron Medium ID (MC)",800,800);
  factEff_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Medium_ele->SetTitle("Factorized efficiency, Electron Medium ID (MC)");
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
  TCanvas *c_factEff_Tight_ele = new TCanvas("c_factEff_Tight_ele","Factorized efficiency, Electron Tight ID (MC)",800,800);
  factEff_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Tight_ele->SetTitle("Factorized efficiency, Electron Tight ID (MC)");
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
  eff_BaselineToLoose_electronTnP  ->Write("eff_Veto_ele"); 
  eff_BaselineToMedium_electronTnP ->Write("eff_Loose_ele"); 
  eff_BaselineToTight_electronTnP  ->Write("eff_Medium_ele");  
  eff_BaselineToVeto_electronTnP   ->Write("eff_Tight_ele");  
  factEff_Veto_ele   ->Write("factEff_Veto_ele");
  factEff_Loose_ele  ->Write("factEff_Loose_ele");
  factEff_Medium_ele ->Write("factEff_Medium_ele");
  factEff_Tight_ele  ->Write("factEff_Tight_ele"); 
  efficiencies_file->Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void factorize_data_ele(string plots_dir, string root_dir) {
  Float_t ele_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t ele_eta_bins[] = {0, 1.479, 2.4};
  Int_t n_ele_pt_bins=12;
  Int_t n_ele_eta_bins=2;
  Float_t mu_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t mu_eta_bins[] = {0, 1.5, 2.1, 2.4};
  Int_t n_mu_pt_bins=12;
  Int_t n_mu_eta_bins=3;
  
  TFile *f_BaselineToLoose_electronTnP          = TFile::Open((plots_dir+"SingleElectron+Run2015D_BaselineToLoose_electronTnP/eff.root"       ).c_str() ,"READ");
  TFile *f_BaselineToMedium_electronTnP         = TFile::Open((plots_dir+"SingleElectron+Run2015D_BaselineToMedium_electronTnP/eff.root"      ).c_str() ,"READ");
  TFile *f_BaselineToTight_electronTnP          = TFile::Open((plots_dir+"SingleElectron+Run2015D_BaselineToTight_electronTnP/eff.root"       ).c_str() ,"READ");
  TFile *f_BaselineToVeto_electronTnP           = TFile::Open((plots_dir+"SingleElectron+Run2015D_BaselineToVeto_electronTnP/eff.root"        ).c_str() ,"READ");
  TFile *f_LooseIdToLooseIdIso_electronTnP      = TFile::Open((plots_dir+"SingleElectron+Run2015D_LooseIdToLooseIdIso_electronTnP/eff.root"   ).c_str() ,"READ");
  TFile *f_LooseIsoToLooseIdIso_electronTnP     = TFile::Open((plots_dir+"SingleElectron+Run2015D_LooseIsoToLooseIdIso_electronTnP/eff.root"  ).c_str() ,"READ");
  TFile *f_MediumIdToMediumIdIso_electronTnP    = TFile::Open((plots_dir+"SingleElectron+Run2015D_MediumIdToMediumIdIso_electronTnP/eff.root" ).c_str() ,"READ");
  TFile *f_MediumIsoToMediumIdIso_electronTnP   = TFile::Open((plots_dir+"SingleElectron+Run2015D_MediumIsoToMediumIdIso_electronTnP/eff.root").c_str() ,"READ");
  TFile *f_TightIdToTightIdIso_electronTnP      = TFile::Open((plots_dir+"SingleElectron+Run2015D_TightIdToTightIdIso_electronTnP/eff.root"   ).c_str() ,"READ");
  TFile *f_TightIsoToTightIdIso_electronTnP     = TFile::Open((plots_dir+"SingleElectron+Run2015D_TightIsoToTightIdIso_electronTnP/eff.root"  ).c_str() ,"READ");
  TFile *f_VetoIdToVetoIdIso_electronTnP        = TFile::Open((plots_dir+"SingleElectron+Run2015D_VetoIdToVetoIdIso_electronTnP/eff.root"     ).c_str() ,"READ");
  TFile *f_VetoIsoToVetoIdIso_electronTnP       = TFile::Open((plots_dir+"SingleElectron+Run2015D_VetoIsoToVetoIdIso_electronTnP/eff.root"    ).c_str() ,"READ");
  
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
  TCanvas *c_absdiff_Veto_ele = new TCanvas("c_absdiff_Veto_ele","Difference from factorized efficiency, Electron Veto ID (Data)",800,800);
  absdiff_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Veto_ele->SetTitle("Difference from factorized efficiency, Electron Veto ID (Data)");
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
  c_absdiff_Veto_ele->Print((plots_dir+"SingleElectron+Run2015D_absdiff_Veto_ele.png").c_str());
  TCanvas *c_absdiff_Loose_ele = new TCanvas("c_absdiff_Loose_ele","Difference from factorized efficiency, Electron Loose ID (Data)",800,800);
  absdiff_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Loose_ele->SetTitle("Difference from factorized efficiency, Electron Loose ID (Data)");
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
  c_absdiff_Loose_ele->Print((plots_dir+"SingleElectron+Run2015D_absdiff_Loose_ele.png").c_str());
  TCanvas *c_absdiff_Medium_ele = new TCanvas("c_absdiff_Medium_ele","Difference from factorized efficiency, Electron Medium ID (Data)",800,800);
  absdiff_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Medium_ele->SetTitle("Difference from factorized efficiency, Electron Medium ID (Data)");
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
  c_absdiff_Medium_ele->Print((plots_dir+"SingleElectron+Run2015D_absdiff_Medium_ele.png").c_str());
  TCanvas *c_absdiff_Tight_ele = new TCanvas("c_absdiff_Tight_ele","Difference from factorized efficiency, Electron Tight ID (Data)",800,800);
  absdiff_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absdiff_Tight_ele->SetTitle("Difference from factorized efficiency, Electron Tight ID (Data)");
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
  c_absdiff_Tight_ele->Print((plots_dir+"SingleElectron+Run2015D_absdiff_Tight_ele.png").c_str());

  // Draw unfactorized efficiency plots
  TCanvas *c_eff_Veto_ele = new TCanvas("c_eff_Veto_ele","Unfactorized efficiency, Electron Veto ID (Data)",800,800);
  eff_BaselineToVeto_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToVeto_electronTnP->SetTitle("Unfactorized efficiency, Electron Veto ID (Data)");
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
  c_eff_Veto_ele->Print((plots_dir+"SingleElectron+Run2015D_eff_Veto_ele.png").c_str());
  TCanvas *c_eff_Loose_ele = new TCanvas("c_eff_Loose_ele","Unfactorized efficiency, Electron Loose ID (Data)",800,800);
  eff_BaselineToLoose_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToLoose_electronTnP->SetTitle("Unfactorized efficiency, Electron Loose ID (Data)");
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
  c_eff_Loose_ele->Print((plots_dir+"SingleElectron+Run2015D_eff_Loose_ele.png").c_str());
  TCanvas *c_eff_Medium_ele = new TCanvas("c_eff_Medium_ele","Unfactorized efficiency, Electron Medium ID (Data)",800,800);
  eff_BaselineToMedium_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToMedium_electronTnP->SetTitle("Unfactorized efficiency, Electron Medium ID (Data)");
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
  c_eff_Medium_ele->Print((plots_dir+"SingleElectron+Run2015D_eff_Medium_ele.png").c_str());
  TCanvas *c_eff_Tight_ele = new TCanvas("c_eff_Tight_ele","Unfactorized efficiency, Electron Tight ID (Data)",800,800);
  eff_BaselineToTight_electronTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToTight_electronTnP->SetTitle("Unfactorized efficiency, Electron Tight ID (Data)");
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
  c_eff_Tight_ele->Print((plots_dir+"SingleElectron+Run2015D_eff_Tight_ele.png").c_str());

  // Draw factorized efficiency plots
  TCanvas *c_factEff_Veto_ele = new TCanvas("c_factEff_Veto_ele","Factorized efficiency, Electron Veto ID (Data)",800,800);
  factEff_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Veto_ele->SetTitle("Factorized efficiency, Electron Veto ID (Data)");
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
  c_factEff_Veto_ele->Print((plots_dir+"SingleElectron+Run2015D_factEff_Veto_ele.png").c_str());
  TCanvas *c_factEff_Loose_ele = new TCanvas("c_factEff_Loose_ele","Factorized efficiency, Electron Loose ID (Data)",800,800);
  factEff_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Loose_ele->SetTitle("Factorized efficiency, Electron Loose ID (Data)");
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
  c_factEff_Loose_ele->Print((plots_dir+"SingleElectron+Run2015D_factEff_Loose_ele.png").c_str());
  TCanvas *c_factEff_Medium_ele = new TCanvas("c_factEff_Medium_ele","Factorized efficiency, Electron Medium ID (Data)",800,800);
  factEff_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Medium_ele->SetTitle("Factorized efficiency, Electron Medium ID (Data)");
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
  c_factEff_Medium_ele->Print((plots_dir+"SingleElectron+Run2015D_factEff_Medium_ele.png").c_str());
  TCanvas *c_factEff_Tight_ele = new TCanvas("c_factEff_Tight_ele","Factorized efficiency, Electron Tight ID (Data)",800,800);
  factEff_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factEff_Tight_ele->SetTitle("Factorized efficiency, Electron Tight ID (Data)");
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
  c_factEff_Tight_ele->Print((plots_dir+"SingleElectron+Run2015D_factEff_Tight_ele.png").c_str());
  
  TFile *efficiencies_file = TFile::Open((root_dir+"SingleElectron+Run2015D_efficiencies_electronTnP.root").c_str(),"RECREATE");
  eff_BaselineToLoose_electronTnP  ->Write("eff_Veto_ele"); 
  eff_BaselineToMedium_electronTnP ->Write("eff_Loose_ele"); 
  eff_BaselineToTight_electronTnP  ->Write("eff_Medium_ele");  
  eff_BaselineToVeto_electronTnP   ->Write("eff_Tight_ele");  
  factEff_Veto_ele   ->Write("factEff_Veto_ele");
  factEff_Loose_ele  ->Write("factEff_Loose_ele");
  factEff_Medium_ele ->Write("factEff_Medium_ele");
  factEff_Tight_ele  ->Write("factEff_Tight_ele"); 
  efficiencies_file->Close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void efficiencies_mc_mu(string plots_dir, string root_dir) {
  Float_t ele_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t ele_eta_bins[] = {0, 1.479, 2.4};
  Int_t n_ele_pt_bins=12;
  Int_t n_ele_eta_bins=2;
  Float_t mu_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t mu_eta_bins[] = {0, 1.5, 2.1, 2.4};
  Int_t n_mu_pt_bins=12;
  Int_t n_mu_eta_bins=3;
  
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
  TCanvas *c_eff_Veto_mu = new TCanvas("c_eff_Veto_mu","Unfactorized efficiency, Muon Veto ID (MC)",800,800);
  eff_BaselineToVeto_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToVeto_muonTnP->SetTitle("Unfactorized efficiency, Muon Veto ID (MC)");
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
  TCanvas *c_eff_Loose_mu = new TCanvas("c_eff_Loose_mu","Unfactorized efficiency, Muon Loose ID (MC)",800,800);
  eff_BaselineToLoose_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToLoose_muonTnP->SetTitle("Unfactorized efficiency, Muon Loose ID (MC)");
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
  TCanvas *c_eff_Medium_mu = new TCanvas("c_eff_Medium_mu","Unfactorized efficiency, Muon Medium ID (MC)",800,800);
  eff_BaselineToMedium_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToMedium_muonTnP->SetTitle("Unfactorized efficiency, Muon Medium ID (MC)");
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
  TCanvas *c_eff_Tight_mu = new TCanvas("c_eff_Tight_mu","Unfactorized efficiency, Muon Tight ID (MC)",800,800);
  eff_BaselineToTight_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToTight_muonTnP->SetTitle("Unfactorized efficiency, Muon Tight ID (MC)");
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
  eff_BaselineToLoose_muonTnP  ->Write("eff_Veto_mu"); 
  eff_BaselineToMedium_muonTnP ->Write("eff_Loose_mu"); 
  eff_BaselineToTight_muonTnP  ->Write("eff_Medium_mu");  
  eff_BaselineToVeto_muonTnP   ->Write("eff_Tight_mu");  
  efficiencies_file->Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void efficiencies_data_mu(string plots_dir, string root_dir) {
  Float_t ele_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t ele_eta_bins[] = {0, 1.479, 2.4};
  Int_t n_ele_pt_bins=12;
  Int_t n_ele_eta_bins=2;
  Float_t mu_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t mu_eta_bins[] = {0, 1.5, 2.1, 2.4};
  Int_t n_mu_pt_bins=12;
  Int_t n_mu_eta_bins=3;
  
  TFile *f_BaselineToLoose_muonTnP          = TFile::Open((plots_dir+"SingleMuon+Run2015D_BaselineToLoose_muonTnP/eff.root" ).c_str()       ,"READ");
  TFile *f_BaselineToMedium_muonTnP         = TFile::Open((plots_dir+"SingleMuon+Run2015D_BaselineToMedium_muonTnP/eff.root").c_str()       ,"READ");
  TFile *f_BaselineToTight_muonTnP          = TFile::Open((plots_dir+"SingleMuon+Run2015D_BaselineToTight_muonTnP/eff.root" ).c_str()       ,"READ");
  TFile *f_BaselineToVeto_muonTnP           = TFile::Open((plots_dir+"SingleMuon+Run2015D_BaselineToVeto_muonTnP/eff.root"  ).c_str()       ,"READ");
  
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
  TCanvas *c_eff_Veto_mu = new TCanvas("c_eff_Veto_mu","Unfactorized efficiency, Muon Veto ID (Data)",800,800);
  eff_BaselineToVeto_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToVeto_muonTnP->SetTitle("Unfactorized efficiency, Muon Veto ID (Data)");
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
  c_eff_Veto_mu->Print((plots_dir+"SingleMuon+Run2015D_eff_Veto_mu.png").c_str());
  TCanvas *c_eff_Loose_mu = new TCanvas("c_eff_Loose_mu","Unfactorized efficiency, Muon Loose ID (Data)",800,800);
  eff_BaselineToLoose_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToLoose_muonTnP->SetTitle("Unfactorized efficiency, Muon Loose ID (Data)");
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
  c_eff_Loose_mu->Print((plots_dir+"SingleMuon+Run2015D_eff_Loose_mu.png").c_str());
  TCanvas *c_eff_Medium_mu = new TCanvas("c_eff_Medium_mu","Unfactorized efficiency, Muon Medium ID (Data)",800,800);
  eff_BaselineToMedium_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToMedium_muonTnP->SetTitle("Unfactorized efficiency, Muon Medium ID (Data)");
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
  c_eff_Medium_mu->Print((plots_dir+"SingleMuon+Run2015D_eff_Medium_mu.png").c_str());
  TCanvas *c_eff_Tight_mu = new TCanvas("c_eff_Tight_mu","Unfactorized efficiency, Muon Tight ID (Data)",800,800);
  eff_BaselineToTight_muonTnP->Draw("TEXTE COLZ");
  gPad->Update(); 
  eff_BaselineToTight_muonTnP->SetTitle("Unfactorized efficiency, Muon Tight ID (Data)");
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
  c_eff_Tight_mu->Print((plots_dir+"SingleMuon+Run2015D_eff_Tight_mu.png").c_str());

  TFile *efficiencies_file = TFile::Open((root_dir+"SingleMuon+Run2015D_efficiencies_muonTnP.root").c_str(),"RECREATE");
  eff_BaselineToLoose_muonTnP  ->Write("eff_Veto_mu"); 
  eff_BaselineToMedium_muonTnP ->Write("eff_Loose_mu"); 
  eff_BaselineToTight_muonTnP  ->Write("eff_Medium_mu");  
  eff_BaselineToVeto_muonTnP   ->Write("eff_Tight_mu");  
  efficiencies_file->Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void make_ele_scale_factors(string plots_dir, string root_dir) {
  Float_t ele_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t ele_eta_bins[] = {0, 1.479, 2.4};
  Int_t n_ele_pt_bins=12;
  Int_t n_ele_eta_bins=2;

  TFile *f_mc   = TFile::Open((root_dir+"DYJetsToLL_efficiencies_electronTnP.root").c_str(),"READ");
  TFile *f_data = TFile::Open((root_dir+"SingleElectron+Run2015D_efficiencies_electronTnP.root").c_str(),"READ");

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

  // Draw unfactorized scalefactor plots
  TCanvas *c_unfactorized_scalefactors_Veto_ele = new TCanvas("c_unfactorized_scalefactors_Veto_ele","Unfactorized scale factors, Electron Veto ID (Data/MC)",800,800);
  unfactorized_scalefactors_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Veto_ele->SetTitle("Unfactorized scale factors, Electron Veto ID (Data/MC)");
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
  TCanvas *c_unfactorized_scalefactors_Loose_ele = new TCanvas("c_unfactorized_scalefactors_Loose_ele","Unfactorized scale factors, Electron Loose ID (Data/MC)",800,800);
  unfactorized_scalefactors_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Loose_ele->SetTitle("Unfactorized scale factors, Electron Loose ID (Data/MC)");
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
  TCanvas *c_unfactorized_scalefactors_Medium_ele = new TCanvas("c_unfactorized_scalefactors_Medium_ele","Unfactorized scale factors, Electron Medium ID (Data/MC)",800,800);
  unfactorized_scalefactors_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Medium_ele->SetTitle("Unfactorized scale factors, Electron Medium ID (Data/MC)");
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
  TCanvas *c_unfactorized_scalefactors_Tight_ele = new TCanvas("c_unfactorized_scalefactors_Tight_ele","Unfactorized scale factors, Electron Tight ID (Data/MC)",800,800);
  unfactorized_scalefactors_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Tight_ele->SetTitle("Unfactorized scale factors, Electron Tight ID (Data/MC)");
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
  TCanvas *c_factorized_scalefactors_Veto_ele = new TCanvas("c_factorized_scalefactors_Veto_ele","Factorized scale factors, Electron Veto ID (Data/MC)",800,800);
  factorized_scalefactors_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factorized_scalefactors_Veto_ele->SetTitle("Factorized scale factors, Electron Veto ID (Data/MC)");
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
  TCanvas *c_factorized_scalefactors_Loose_ele = new TCanvas("c_factorized_scalefactors_Loose_ele","Factorized scale factors, Electron Loose ID (Data/MC)",800,800);
  factorized_scalefactors_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factorized_scalefactors_Loose_ele->SetTitle("Factorized scale factors, Electron Loose ID (Data/MC)");
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
  TCanvas *c_factorized_scalefactors_Medium_ele = new TCanvas("c_factorized_scalefactors_Medium_ele","Factorized scale factors, Electron Medium ID (Data/MC)",800,800);
  factorized_scalefactors_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factorized_scalefactors_Medium_ele->SetTitle("Factorized scale factors, Electron Medium ID (Data/MC)");
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
  TCanvas *c_factorized_scalefactors_Tight_ele = new TCanvas("c_factorized_scalefactors_Tight_ele","Factorized scale factors, Electron Tight ID (Data/MC)",800,800);
  factorized_scalefactors_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  factorized_scalefactors_Tight_ele->SetTitle("Factorized scale factors, Electron Tight ID (Data/MC)");
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

  TFile *scalefactors_file = TFile::Open((root_dir+"scalefactors_ele.root").c_str(),"RECREATE");
  unfactorized_scalefactors_Veto_ele   ->Write("unfactorized_scalefactors_Veto_ele"   );
  unfactorized_scalefactors_Loose_ele  ->Write("unfactorized_scalefactors_Loose_ele"  );
  unfactorized_scalefactors_Medium_ele ->Write("unfactorized_scalefactors_Medium_ele" );
  unfactorized_scalefactors_Tight_ele  ->Write("unfactorized_scalefactors_Tight_ele"  );
  factorized_scalefactors_Veto_ele    ->Write("factorized_scalefactors_Veto_ele"  );
  factorized_scalefactors_Loose_ele   ->Write("factorized_scalefactors_Loose_ele" );
  factorized_scalefactors_Medium_ele  ->Write("factorized_scalefactors_Medium_ele");
  factorized_scalefactors_Tight_ele   ->Write("factorized_scalefactors_Tight_ele" );
  scalefactors_file->Close();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void make_mu_scale_factors(string plots_dir, string root_dir) {
  Float_t mu_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t mu_eta_bins[] = {0, 1.5, 2.1, 2.4};
  Int_t n_mu_pt_bins=12;
  Int_t n_mu_eta_bins=3;
  
  TFile *f_mc   = TFile::Open((root_dir+"DYJetsToLL_efficiencies_muonTnP.root").c_str(),"READ");
  TFile *f_data = TFile::Open((root_dir+"SingleMuon+Run2015D_efficiencies_muonTnP.root").c_str(),"READ");

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

  // Draw unfactorized scalefactor plots
  TCanvas *c_unfactorized_scalefactors_Veto_mu = new TCanvas("c_unfactorized_scalefactors_Veto_mu","Unfactorized scale factors, Muon Veto ID (Data/MC)",800,800);
  unfactorized_scalefactors_Veto_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Veto_mu->SetTitle("Unfactorized scale factors, Muon Veto ID (Data/MC)");
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
  palette_axis = (TPaletteAxis*) unfactorized_scalefactors_Veto_mu->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_unfactorized_scalefactors_Veto_mu->Update();
  c_unfactorized_scalefactors_Veto_mu->Print((plots_dir+"unfactorized_scalefactors_Veto_mu.png").c_str());
  TCanvas *c_unfactorized_scalefactors_Loose_mu = new TCanvas("c_unfactorized_scalefactors_Loose_mu","Unfactorized scale factors, Muon Loose ID (Data/MC)",800,800);
  unfactorized_scalefactors_Loose_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Loose_mu->SetTitle("Unfactorized scale factors, Muon Loose ID (Data/MC)");
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
  TCanvas *c_unfactorized_scalefactors_Medium_mu = new TCanvas("c_unfactorized_scalefactors_Medium_mu","Unfactorized scale factors, Muon Medium ID (Data/MC)",800,800);
  unfactorized_scalefactors_Medium_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Medium_mu->SetTitle("Unfactorized scale factors, Muon Medium ID (Data/MC)");
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
  TCanvas *c_unfactorized_scalefactors_Tight_mu = new TCanvas("c_unfactorized_scalefactors_Tight_mu","Unfactorized scale factors, Muon Tight ID (Data/MC)",800,800);
  unfactorized_scalefactors_Tight_mu->Draw("TEXTE COLZ");
  gPad->Update(); 
  unfactorized_scalefactors_Tight_mu->SetTitle("Unfactorized scale factors, Muon Tight ID (Data/MC)");
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

  TFile *scalefactors_file = TFile::Open((root_dir+"scalefactors_mu.root").c_str(),"RECREATE");
  unfactorized_scalefactors_Veto_mu   ->Write("unfactorized_scalefactors_Veto_mu"  );
  unfactorized_scalefactors_Loose_mu  ->Write("unfactorized_scalefactors_Loose_mu" );
  unfactorized_scalefactors_Medium_mu ->Write("unfactorized_scalefactors_Medium_mu");
  unfactorized_scalefactors_Tight_mu  ->Write("unfactorized_scalefactors_Tight_mu" );
  scalefactors_file->Close();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void plot_ele() {
  gStyle->SetOptStat(0);
  Float_t bins[] = {20,30,40,60,80,150};
  Int_t nbins=5;
  
  TH1D *h_data_eff_veto   = new TH1D("h_data_eff_veto","Electron ID efficiencies in 2015C data using tag and probe",nbins,bins);
  TH1D *h_data_eff_medium = new TH1D("h_data_eff_medium","Loose to medium",nbins,bins);
  TH1D *h_data_eff_tight  = new TH1D("h_data_eff_tight","Loose to tight",nbins,bins);
  TH1D *h_mc_eff_veto     = new TH1D("h_mc_eff_veto","Electron ID efficiencies in MC using tag and probe",nbins,bins);
  TH1D *h_mc_eff_medium   = new TH1D("h_mc_eff_medium","Loose to medium",nbins,bins);
  TH1D *h_mc_eff_tight    = new TH1D("h_mc_eff_tight","Loose to tight",nbins,bins);
  TH1D *h_ratio_eff_veto     = new TH1D("h_ratio_eff_veto","Electron ID scale factors from MC to data",nbins,bins);
  TH1D *h_ratio_eff_medium   = new TH1D("h_ratio_eff_medium","Loose to medium",nbins,bins);
  TH1D *h_ratio_eff_tight    = new TH1D("h_ratio_eff_tight","Loose to tight",nbins,bins);

  double data_eff_veto[]       = {0.8272808,  0.9283912,  0.9673420,  0.9699439,  0.9874522};
  double data_eff_medium[]     = {0.6770211,  0.7857058,  0.8492878,  0.8921974,  0.9396030};
  double data_eff_tight[]      = {0.5638646,  0.6700576,  0.7592614,  0.8063168,  0.9065888};
  double data_err_eff_veto[]   = {0.0239439,  0.0052142,  0.0030022,  0.0117742,  0.0133755};
  double data_err_eff_medium[] = {0.0204696,  0.0009449,  0.0050543,  0.0220630,  0.0252899};
  double data_err_eff_tight[]  = {0.0168675,  0.0072819,  0.0055368,  0.0223588,  0.0357888};
  double mc_eff_veto[]         = {0.9410147,  0.9691390,  0.9868144,  0.9906645,  0.9925915};
  double mc_eff_medium[]       = {0.7576817,  0.8388063,  0.8917681,  0.9197389,  0.9359847};
  double mc_eff_tight[]        = {0.6287441,  0.7210381,  0.7966400,  0.8474214,  0.8831588};
  double mc_err_eff_veto[]     = {0.0005299,  0.0001449,  0.0000881,  0.0002378,  0.0002474};
  double mc_err_eff_medium[]   = {0.0006801,  0.0003545,  0.0002198,  0.0005312,  0.0007027};
  double mc_err_eff_tight[]    = {0.0007327,  0.0000890,  0.0002817,  0.0006866,  0.0008913};
  for(int i=1; i<=nbins; i++) {
    h_data_eff_veto   ->SetBinContent (i, data_eff_veto[i-1]       );  
    h_data_eff_veto   ->SetBinError   (i, data_err_eff_veto[i-1]   );  
    h_data_eff_medium ->SetBinContent (i, data_eff_medium[i-1]     );  
    h_data_eff_medium ->SetBinError   (i, data_err_eff_medium[i-1] );  
    h_data_eff_tight  ->SetBinContent (i, data_eff_tight[i-1]      );  
    h_data_eff_tight  ->SetBinError   (i, data_err_eff_tight[i-1]  );  
    h_mc_eff_veto     ->SetBinContent (i, mc_eff_veto[i-1]       );  
    h_mc_eff_veto     ->SetBinError   (i, mc_err_eff_veto[i-1]   );  
    h_mc_eff_medium   ->SetBinContent (i, mc_eff_medium[i-1]     );  
    h_mc_eff_medium   ->SetBinError   (i, mc_err_eff_medium[i-1] );  
    h_mc_eff_tight    ->SetBinContent (i, mc_eff_tight[i-1]      );  
    h_mc_eff_tight    ->SetBinError   (i, mc_err_eff_tight[i-1]  ); 
    // add relative errors in quadrature
    h_ratio_eff_veto   ->SetBinContent(i, data_eff_veto[i-1]/mc_eff_veto[i-1]);
    h_ratio_eff_veto   ->SetBinError  (i, (data_eff_veto[i-1]/mc_eff_veto[i-1]) * sqrt(pow(data_err_eff_veto[i-1]/data_eff_veto[i-1], 2) + pow(mc_err_eff_veto[i-1]/mc_eff_veto[i-1],2) ));
    h_ratio_eff_medium ->SetBinContent(i, data_eff_medium[i-1]/mc_eff_medium[i-1]);
    h_ratio_eff_medium ->SetBinError  (i, (data_eff_medium[i-1]/mc_eff_medium[i-1]) * sqrt(pow(data_err_eff_medium[i-1]/data_eff_medium[i-1], 2) + pow(mc_err_eff_medium[i-1]/mc_eff_medium[i-1],2) ));
    h_ratio_eff_tight   ->SetBinContent(i, data_eff_tight[i-1]/mc_eff_tight[i-1]);
    h_ratio_eff_tight   ->SetBinError  (i, (data_eff_tight[i-1]/mc_eff_tight[i-1]) * sqrt(pow(data_err_eff_tight[i-1]/data_eff_tight[i-1], 2) + pow(mc_err_eff_tight[i-1]/mc_eff_tight[i-1],2) ));
    printf("pT [ %d , %d ] GeV: scale factor %f +/- %f\n", int(bins[i-1]), int(bins[i]), h_ratio_eff_tight->GetBinContent(i), h_ratio_eff_tight->GetBinError(i));
  }
  TCanvas *c1 = new TCanvas("c1","Electron efficiency in data vs. pT");
  h_data_eff_veto   ->SetMarkerStyle(20);
  h_data_eff_medium ->SetMarkerStyle(21);
  h_data_eff_tight  ->SetMarkerStyle(22);
  h_data_eff_veto   ->SetLineColor(52);
  h_data_eff_medium ->SetLineColor(62);
  h_data_eff_tight  ->SetLineColor(98);
  h_data_eff_veto   ->SetMarkerColor(52);
  h_data_eff_medium ->SetMarkerColor(62);
  h_data_eff_tight  ->SetMarkerColor(98);
  h_data_eff_veto->Draw("E1 P");
  h_data_eff_medium->Draw("E1 P SAME");
  h_data_eff_tight->Draw("E1 P SAME");
  h_data_eff_veto->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_data_eff_veto->GetYaxis()->SetTitle("Efficiency");
  h_data_eff_veto->SetMaximum(1.2);
  h_data_eff_veto->SetMinimum(0);
  TLegend *data_legend = new TLegend(0.65,0.15,0.85,0.35);
  data_legend->AddEntry(h_data_eff_veto,"Loose to veto","lp");
  data_legend->AddEntry(h_data_eff_medium,"Loose to medium","lp");
  data_legend->AddEntry(h_data_eff_tight,"Loose to tight","lp");
  data_legend->SetFillColor(0);
  data_legend->Draw("SAME");

  TCanvas *c2 = new TCanvas("c2","Electron efficiency in MC vs. pT");
  h_mc_eff_veto   ->SetMarkerStyle(20);
  h_mc_eff_medium ->SetMarkerStyle(21);
  h_mc_eff_tight  ->SetMarkerStyle(22);
  h_mc_eff_veto   ->SetLineColor(52);
  h_mc_eff_medium ->SetLineColor(62);
  h_mc_eff_tight  ->SetLineColor(98);
  h_mc_eff_veto   ->SetMarkerColor(52);
  h_mc_eff_medium ->SetMarkerColor(62);
  h_mc_eff_tight  ->SetMarkerColor(98);
  h_mc_eff_veto->Draw("E1 P");
  h_mc_eff_medium->Draw("E1 P SAME");
  h_mc_eff_tight->Draw("E1 P SAME");
  h_mc_eff_veto->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_mc_eff_veto->GetYaxis()->SetTitle("Efficiency");
  h_mc_eff_veto->SetMaximum(1.2);
  h_mc_eff_veto->SetMinimum(0);
  TLegend *mc_legend = new TLegend(0.65,0.15,0.85,0.35);
  mc_legend->AddEntry(h_data_eff_veto,"Loose to veto","lp");
  mc_legend->AddEntry(h_data_eff_medium,"Loose to medium","lp");
  mc_legend->AddEntry(h_data_eff_tight,"Loose to tight","lp");
  mc_legend->SetFillColor(0);
  mc_legend->Draw("SAME");

  TCanvas *c3 = new TCanvas("c3","Ratio of efficiency data:MC vs. pT");
  h_ratio_eff_veto   ->SetMarkerStyle(20);
  h_ratio_eff_medium ->SetMarkerStyle(21);
  h_ratio_eff_tight  ->SetMarkerStyle(22);
  h_ratio_eff_veto   ->SetLineColor(52);
  h_ratio_eff_medium ->SetLineColor(62);
  h_ratio_eff_tight  ->SetLineColor(98);
  h_ratio_eff_veto   ->SetMarkerColor(52);
  h_ratio_eff_medium ->SetMarkerColor(62);
  h_ratio_eff_tight  ->SetMarkerColor(98);
  h_ratio_eff_veto->Draw("E1 P");
  h_ratio_eff_medium->Draw("E1 P SAME");
  h_ratio_eff_tight->Draw("E1 P SAME");
  h_ratio_eff_veto->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_ratio_eff_veto->GetYaxis()->SetTitle("Efficiency");
  h_ratio_eff_veto->SetMaximum(1.2);
  h_ratio_eff_veto->SetMinimum(0.8);
  TLegend *ratio_legend = new TLegend(0.65,0.15,0.85,0.35);
  ratio_legend->AddEntry(h_data_eff_veto,"Loose to veto","lp");
  ratio_legend->AddEntry(h_data_eff_medium,"Loose to medium","lp");
  ratio_legend->AddEntry(h_data_eff_tight,"Loose to tight","lp");
  ratio_legend->SetFillColor(0);
  ratio_legend->Draw("SAME");
}
void plot_mu() {
  gStyle->SetOptStat(0);
  Float_t bins[] = {20,30,40,60,80,150};
  Int_t nbins=5;
  
  TH1D *h_eff_tight   = new TH1D("h_eff_tight","Muon ID efficiency using tag and probe",nbins,bins);

  double eff_tight[]  = {0.6110737,  0.7796324,  0.8833061,  0.9420620,  0.9509638};
  double err_eff_tight[]  = {0.0078726,  0.0044903,  0.0030520,  0.0077855,  0.0107436 };
  for(int i=1; i<=nbins; i++) {
    h_eff_tight->SetBinContent(i, eff_tight[i-1]);  
    h_eff_tight->SetBinError(i, err_eff_tight[i-1]);  
  }
  TCanvas *c1 = new TCanvas("c1","Muon efficiency as function of pT");
  //h_eff_veto   ->SetMarkerStyle(20);
  //h_eff_medium ->SetMarkerStyle(21);
  h_eff_tight  ->SetMarkerStyle(22);
  //h_eff_veto   ->SetLineColor(52);
  //h_eff_medium ->SetLineColor(62);
  h_eff_tight  ->SetLineColor(98);
  //h_eff_veto   ->SetMarkerColor(52);
  //h_eff_medium ->SetMarkerColor(62);
  h_eff_tight  ->SetMarkerColor(98);
  h_eff_tight->Draw("E1 P");
  //h_eff_veto->Draw("E1 P");
  //h_eff_medium->Draw("E1 P SAME");
  //h_eff_tight->Draw("E1 P SAME");
  h_eff_tight->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_eff_tight->GetYaxis()->SetTitle("Efficiency");
  h_eff_tight->SetMaximum(1.2);
  h_eff_tight->SetMinimum(0);
  TLegend *legend = new TLegend(0.65,0.15,0.85,0.35);
  //legend->AddEntry(h_eff_veto,"Loose to veto","lp");
  //legend->AddEntry(h_eff_medium,"Loose to medium","lp");
  legend->AddEntry(h_eff_tight,"Loose to tight","lp");
  legend->SetFillColor(0);
  legend->Draw("SAME");
}



