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
#include <TLegend.h>
#include <TFile.h>
#include <TCut.h>
#include <leptons.h>
#include <TColor.h>
#include <TPaletteAxis.h>
#include <iostream>
Int_t mit_red  = 1861; 
Int_t mit_gray = 1862; 
Float_t ele_pt_bins[] = {10,20,30,40,50,200};
Float_t ele_eta_bins[] = {0., 0.8, 1.4442, 1.566, 2.0, 2.5};
Int_t n_ele_pt_bins=5;
Int_t n_ele_eta_bins=5;
Float_t mu_pt_bins[] = {10,20,30,40,50,200,8000};
Float_t mu_eta_bins[] = {0, 0.2, 0.3, 0.9, 1.2, 2.1, 2.4};
Int_t n_mu_pt_bins=6;
Int_t n_mu_eta_bins=6;
int marker_colors[] = {1, mit_red, mit_gray, 9, 6, 4, 8};
int marker_styles[] = {20, 21, 22, 23, 33, 34, 20};
//Float_t ele_pt_bins[] = {10,20,30,40,50,70,100,8000};
//Float_t ele_eta_bins[] = {0, 1.479, 2.4};
//Int_t n_ele_pt_bins=7;
//Int_t n_ele_eta_bins=2;
//Float_t mu_pt_bins[] = {10,20,30,40,50,70,100,8000};
//Float_t mu_eta_bins[] = {0, 1.479, 2.4};
//Int_t n_mu_pt_bins=7;
//Int_t n_mu_eta_bins=2;

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
void scale_factors(string plots_dir, string root_dir, string basename_config) {
  // This function calculates the efficiencies based on output from MIT TNP saved in subdirectories of plots_dir
  // The base names of the subdirectories are recorded in the config file
  // The efficiencies are plotted in plots_dir
  // The efficiencies and scale factors are recorded in root_dir in a 
  // rootfile whose filename is taken from basename_config

  // Pad directories with a slash at the end if it's not there
  if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
  if( root_dir[plots_dir.size()-1]  != '/' )  root_dir  = root_dir + "/";

  //unsigned int array_size=32;
  // Set up arrays
  //TObjArray f_mc_            (array_size);
  //TObjArray h_eff_mc_        (array_size);
  //TObjArray h_error_lo_mc_   (array_size);
  //TObjArray h_error_hi_mc_   (array_size);
  //TObjArray c_eff_mc_        (array_size);

  //TObjArray f_data_            (array_size);
  //TObjArray h_eff_data_        (array_size);
  //TObjArray h_error_lo_data_   (array_size);
  //TObjArray h_error_hi_data_   (array_size);
  //TObjArray c_eff_data_        (array_size);

  // Open the output rootfile
  string output_rootfile_name;
  size_t ext_location = basename_config.find_last_of(".");
  if(ext_location != string::npos) output_rootfile_name = "scalefactors_"+basename_config.substr(0, ext_location)+".root";
  else output_rootfile_name = "scalefactors_"+basename_config+".root";
  TFile *output_rootfile = TFile::Open(output_rootfile_name.c_str(), "RECREATE");

  // Now read config file. Each line has
  // <data basename> <mc basename> <selection> <flavor>
  ifstream config_stream;
  config_stream.open(basename_config.c_str());
  assert(config_stream.is_open());
  string line, data_basename, mc_basename, output_basename, selection, flavor;
  while(getline(config_stream, line)) {
    stringstream ss(line);
    ss >> data_basename >> mc_basename >> selection >> flavor;

    if(flavor == "electron")  output_basename = selection+"_ele";
    else if(flavor == "muon")   output_basename = selection+"_mu";
    else                    output_basename = selection+"_"+flavor;
    
    TFile *f_data = TFile::Open( ( plots_dir + data_basename   + "/eff.root").c_str(), "READ");
    TH2D *h_eff_data = (TH2D*) f_data->Get("hEffEtaPt");
    TH2D *h_error_lo_data = (TH2D*) f_data->Get("hErrlEtaPt");
    TH2D *h_error_hi_data = (TH2D*) f_data->Get("hErrhEtaPt");
    h_eff_data->SetName(("eff_data_"+output_basename).c_str());
    h_eff_data->SetTitle(("Efficiency for "+selection+" "+flavor+" selection (Data)").c_str());
    
    TFile *f_mc   = TFile::Open( ( plots_dir + mc_basename + "/eff.root").c_str(), "READ");
    TH2D *h_eff_mc   = (TH2D*) f_mc->Get("hEffEtaPt");
    TH2D *h_error_lo_mc = (TH2D*) f_mc->Get("hErrlEtaPt");
    TH2D *h_error_hi_mc = (TH2D*) f_mc->Get("hErrhEtaPt");
    h_eff_mc->SetName(("eff_mc_"+output_basename).c_str());
    h_eff_mc->SetTitle(("Efficiency for "+selection+" "+flavor+" selection (MC)").c_str());
    
    // Divide Data/MC to get the scale factors
    TH2D *h_sf = (TH2D*) h_eff_data->Clone();
    h_sf->SetName(("scalefactors_"+output_basename).c_str());
    h_sf->SetTitle(("Scale factors for "+selection+" "+flavor+" selection").c_str());

    // Propagate asymmetric statistical errors
    TH2D *h_sf_error_lo = (TH2D*) h_sf->Clone();
    TH2D *h_sf_error_hi = (TH2D*) h_sf->Clone();
    h_sf_error_lo->SetName(("scalefactors_"+output_basename+"_stat_error_lo").c_str());
    h_sf_error_hi->SetName(("scalefactors_"+output_basename+"_stat_error_hi").c_str());
    h_sf_error_lo->SetTitle(("scalefactors_"+output_basename+"_stat_error_lo").c_str());
    h_sf_error_hi->SetTitle(("scalefactors_"+output_basename+"_stat_error_hi").c_str());

    for(int i = 1; i <= h_eff_mc->GetNbinsX(); i++) { for(int j = 1; j <= h_eff_mc->GetNbinsY(); j++) {
      unsigned int nbin = h_eff_mc->GetBin(i,j);
      h_sf->SetBinContent(nbin, h_eff_data->GetBinContent(nbin) / h_eff_mc->GetBinContent(nbin));
      h_sf_error_hi->SetBinContent(nbin,
        h_sf->GetBinContent(nbin) * sqrt(
          pow(h_error_hi_data->GetBinContent(nbin) / h_eff_data->GetBinContent(nbin), 2) +
          pow(h_error_lo_mc->GetBinContent(nbin) / h_eff_mc->GetBinContent(nbin),2)
        )
      );
      h_sf_error_lo->SetBinContent(nbin,
        h_sf->GetBinContent(nbin) * sqrt(
          pow(h_error_lo_data->GetBinContent(nbin) / h_eff_data->GetBinContent(nbin), 2) +
          pow(h_error_hi_mc->GetBinContent(nbin) / h_eff_mc->GetBinContent(nbin),2)
        )
      );
      // Choose the max for the s.f. histogram errors
      // Analyzers who naively use this uncertainty value will just get worse sensitivity :^)
      h_sf->SetBinError(nbin, TMath::Max(
        h_sf_error_hi->GetBinContent(nbin),
        h_sf_error_lo->GetBinContent(nbin)
      ));
    }} 

    // Start drawing stuff 
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.3f");
    mitPalette();
    TPaletteAxis *palette_axis;

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800,600);
    h_eff_data->Draw("TEXTE COLZ");
    canvas->Update();
    h_eff_data->GetXaxis()->SetTitle("| #eta |");
    h_eff_data->GetXaxis()->SetTitleOffset(0.9);
    h_eff_data->GetXaxis()->SetTitleSize(0.04);
    h_eff_data->GetXaxis()->SetLabelSize(0.02);
    h_eff_data->GetYaxis()->SetTitle("p_{T} [GeV]");
    h_eff_data->GetYaxis()->SetTitleOffset(0.9);
    h_eff_data->GetYaxis()->SetTitleSize(0.04);
    h_eff_data->GetYaxis()->SetLabelSize(0.02);
    h_eff_data->GetYaxis()->SetRangeUser(10,200);
    h_eff_data->SetMinimum(0);
    h_eff_data->SetMaximum(1.);
    h_eff_data->SetMarkerSize(.9);
    palette_axis = (TPaletteAxis*) h_eff_data->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    canvas->Update();
    canvas->Print((plots_dir + string(h_eff_data->GetName()) + ".png").c_str());

    h_eff_mc->Draw("TEXTE COLZ");
    canvas->Update();
    h_eff_mc->GetXaxis()->SetTitle("| #eta |");
    h_eff_mc->GetXaxis()->SetTitleOffset(0.9);
    h_eff_mc->GetXaxis()->SetTitleSize(0.04);
    h_eff_mc->GetXaxis()->SetLabelSize(0.02);
    h_eff_mc->GetYaxis()->SetTitle("p_{T} [GeV]");
    h_eff_mc->GetYaxis()->SetTitleOffset(0.9);
    h_eff_mc->GetYaxis()->SetTitleSize(0.04);
    h_eff_mc->GetYaxis()->SetLabelSize(0.02);
    h_eff_mc->GetYaxis()->SetRangeUser(10,200);
    h_eff_mc->SetMinimum(0);
    h_eff_mc->SetMaximum(1.);
    h_eff_mc->SetMarkerSize(.9);
    palette_axis = (TPaletteAxis*) h_eff_mc->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    canvas->Update();
    canvas->Print((plots_dir + string(h_eff_mc->GetName()) + ".png").c_str());

    h_sf->Draw("TEXTE COLZ");
    canvas->Update();
    h_sf->GetXaxis()->SetTitle("| #eta |");
    h_sf->GetXaxis()->SetTitleOffset(0.9);
    h_sf->GetXaxis()->SetTitleSize(0.04);
    h_sf->GetXaxis()->SetLabelSize(0.02);
    h_sf->GetYaxis()->SetTitle("p_{T} [GeV]");
    h_sf->GetYaxis()->SetTitleOffset(0.9);
    h_sf->GetYaxis()->SetTitleSize(0.04);
    h_sf->GetYaxis()->SetLabelSize(0.02);
    h_sf->GetYaxis()->SetRangeUser(10,200);
    h_sf->SetMinimum(.8);
    h_sf->SetMaximum(1.2);
    h_sf->SetMarkerSize(.9);
    palette_axis = (TPaletteAxis*) h_sf->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    canvas->Update();
    canvas->Print((plots_dir + string(h_sf->GetName()) + ".png").c_str());

    // Write to file
    output_rootfile->cd();
    h_eff_data->Write();
    h_eff_mc->Write();
    h_sf->Write();
    h_sf_error_lo->Write();
    h_sf_error_hi->Write();
    
    delete canvas;
    delete h_eff_data;
    delete h_eff_mc;
    delete h_sf;
    delete h_sf_error_lo;
    delete h_sf_error_hi;
    f_data->Close();
    f_mc->Close();
  }
  output_rootfile->Close();
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
  //TH2D *factorized_scalefactors_Veto_ele_method1   = (TH2D*)ele_scalefactors_file_method1->Get("factorized_scalefactors_Veto_ele");
  //TH2D *factorized_scalefactors_Loose_ele_method1  = (TH2D*)ele_scalefactors_file_method1->Get("factorized_scalefactors_Loose_ele");
  //TH2D *factorized_scalefactors_Medium_ele_method1 = (TH2D*)ele_scalefactors_file_method1->Get("factorized_scalefactors_Medium_ele");
  //TH2D *factorized_scalefactors_Tight_ele_method1  = (TH2D*)ele_scalefactors_file_method1->Get("factorized_scalefactors_Tight_ele");
  //TH2D *factorized_scalefactors_Veto_ele_method2   = (TH2D*)ele_scalefactors_file_method2->Get("factorized_scalefactors_Veto_ele");
  //TH2D *factorized_scalefactors_Loose_ele_method2  = (TH2D*)ele_scalefactors_file_method2->Get("factorized_scalefactors_Loose_ele");
  //TH2D *factorized_scalefactors_Medium_ele_method2 = (TH2D*)ele_scalefactors_file_method2->Get("factorized_scalefactors_Medium_ele");
  //TH2D *factorized_scalefactors_Tight_ele_method2  = (TH2D*)ele_scalefactors_file_method2->Get("factorized_scalefactors_Tight_ele");
  TH2D *unfactorized_scalefactors_Veto_ele_method1   = (TH2D*)ele_scalefactors_file_method1->Get("unfactorized_scalefactors_Veto_ele");
  TH2D *unfactorized_scalefactors_Loose_ele_method1  = (TH2D*)ele_scalefactors_file_method1->Get("unfactorized_scalefactors_Loose_ele");
  TH2D *unfactorized_scalefactors_Medium_ele_method1 = (TH2D*)ele_scalefactors_file_method1->Get("unfactorized_scalefactors_Medium_ele");
  TH2D *unfactorized_scalefactors_Tight_ele_method1  = (TH2D*)ele_scalefactors_file_method1->Get("unfactorized_scalefactors_Tight_ele");
  TH2D *unfactorized_scalefactors_Veto_ele_method2   = (TH2D*)ele_scalefactors_file_method2->Get("unfactorized_scalefactors_Veto_ele");
  TH2D *unfactorized_scalefactors_Loose_ele_method2  = (TH2D*)ele_scalefactors_file_method2->Get("unfactorized_scalefactors_Loose_ele");
  TH2D *unfactorized_scalefactors_Medium_ele_method2 = (TH2D*)ele_scalefactors_file_method2->Get("unfactorized_scalefactors_Medium_ele");
  TH2D *unfactorized_scalefactors_Tight_ele_method2  = (TH2D*)ele_scalefactors_file_method2->Get("unfactorized_scalefactors_Tight_ele");
  
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
  //TH2D *absDiff_factorized_scalefactors_Veto_ele = new TH2D(
  //  "absDiff_factorized_scalefactors_Veto_ele",
  //  "Absolute difference between factorized Electron Veto scalefactors computed 2 ways",
  //  n_ele_eta_bins, ele_eta_bins,
  //  n_ele_pt_bins, ele_pt_bins
  //);
  //TH2D *absDiff_factorized_scalefactors_Loose_ele = new TH2D(
  //  "absDiff_factorized_scalefactors_Loose_ele",
  //  "Absolute difference between factorized Electron Loose scalefactors computed 2 ways",
  //  n_ele_eta_bins, ele_eta_bins,
  //  n_ele_pt_bins, ele_pt_bins
  //);
  //TH2D *absDiff_factorized_scalefactors_Medium_ele = new TH2D(
  //  "absDiff_factorized_scalefactors_Medium_ele",
  //  "Absolute difference between factorized Electron Medium scalefactors computed 2 ways",
  //  n_ele_eta_bins, ele_eta_bins,
  //  n_ele_pt_bins, ele_pt_bins
  //);
  //TH2D *absDiff_factorized_scalefactors_Tight_ele = new TH2D(
  //  "absDiff_factorized_scalefactors_Tight_ele",
  //  "Absolute difference between factorized Electron Tight scalefactors computed 2 ways",
  //  n_ele_eta_bins, ele_eta_bins,
  //  n_ele_pt_bins, ele_pt_bins
  //);
  TH2D *absDiff_unfactorized_scalefactors_Veto_ele = new TH2D(
    "absDiff_unfactorized_scalefactors_Veto_ele",
    "Absolute difference between unfactorized Electron Veto scalefactors computed 2 ways",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  TH2D *absDiff_unfactorized_scalefactors_Loose_ele = new TH2D(
    "absDiff_unfactorized_scalefactors_Loose_ele",
    "Absolute difference between unfactorized Electron Loose scalefactors computed 2 ways",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  TH2D *absDiff_unfactorized_scalefactors_Medium_ele = new TH2D(
    "absDiff_unfactorized_scalefactors_Medium_ele",
    "Absolute difference between unfactorized Electron Medium scalefactors computed 2 ways",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );
  TH2D *absDiff_unfactorized_scalefactors_Tight_ele = new TH2D(
    "absDiff_unfactorized_scalefactors_Tight_ele",
    "Absolute difference between unfactorized Electron Tight scalefactors computed 2 ways",
    n_ele_eta_bins, ele_eta_bins,
    n_ele_pt_bins, ele_pt_bins
  );

  // compute difference
  (*absDiff_unfactorized_scalefactors_Veto_mu)   = (*unfactorized_scalefactors_Veto_mu_method1)   - (*unfactorized_scalefactors_Veto_mu_method2);
  (*absDiff_unfactorized_scalefactors_Loose_mu)  = (*unfactorized_scalefactors_Loose_mu_method1)  - (*unfactorized_scalefactors_Loose_mu_method2);
  (*absDiff_unfactorized_scalefactors_Medium_mu) = (*unfactorized_scalefactors_Medium_mu_method1) - (*unfactorized_scalefactors_Medium_mu_method2);
  (*absDiff_unfactorized_scalefactors_Tight_mu)  = (*unfactorized_scalefactors_Tight_mu_method1)  - (*unfactorized_scalefactors_Tight_mu_method2);
  //(*absDiff_factorized_scalefactors_Veto_ele)   = (*factorized_scalefactors_Veto_ele_method1)   - (*factorized_scalefactors_Veto_ele_method2);
  //(*absDiff_factorized_scalefactors_Loose_ele)  = (*factorized_scalefactors_Loose_ele_method1)  - (*factorized_scalefactors_Loose_ele_method2);
  //(*absDiff_factorized_scalefactors_Medium_ele) = (*factorized_scalefactors_Medium_ele_method1) - (*factorized_scalefactors_Medium_ele_method2);
  //(*absDiff_factorized_scalefactors_Tight_ele)  = (*factorized_scalefactors_Tight_ele_method1)  - (*factorized_scalefactors_Tight_ele_method2);
  (*absDiff_unfactorized_scalefactors_Veto_ele)   = (*unfactorized_scalefactors_Veto_ele_method1)   - (*unfactorized_scalefactors_Veto_ele_method2);
  (*absDiff_unfactorized_scalefactors_Loose_ele)  = (*unfactorized_scalefactors_Loose_ele_method1)  - (*unfactorized_scalefactors_Loose_ele_method2);
  (*absDiff_unfactorized_scalefactors_Medium_ele) = (*unfactorized_scalefactors_Medium_ele_method1) - (*unfactorized_scalefactors_Medium_ele_method2);
  (*absDiff_unfactorized_scalefactors_Tight_ele)  = (*unfactorized_scalefactors_Tight_ele_method1)  - (*unfactorized_scalefactors_Tight_ele_method2);
  // take its absolute value
  for(int i = 1; i<=n_mu_eta_bins; i++) { for(int j = 1; j<=n_mu_pt_bins; j++) {
    int nbin = absDiff_unfactorized_scalefactors_Veto_mu->GetBin(i,j);
    absDiff_unfactorized_scalefactors_Veto_mu   ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Veto_mu  ->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Loose_mu  ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Loose_mu ->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Medium_mu ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Medium_mu->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Tight_mu  ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Tight_mu ->GetBinContent(nbin)));  
  }}
  for(int i = 1; i<=n_ele_eta_bins; i++) { for(int j = 1; j<=n_ele_pt_bins; j++) {
    int nbin = absDiff_unfactorized_scalefactors_Veto_ele->GetBin(i,j);
    //absDiff_factorized_scalefactors_Veto_ele   ->SetBinContent(nbin, TMath::Abs(absDiff_factorized_scalefactors_Veto_ele  ->GetBinContent(nbin)));  
    //absDiff_factorized_scalefactors_Loose_ele  ->SetBinContent(nbin, TMath::Abs(absDiff_factorized_scalefactors_Loose_ele ->GetBinContent(nbin)));  
    //absDiff_factorized_scalefactors_Medium_ele ->SetBinContent(nbin, TMath::Abs(absDiff_factorized_scalefactors_Medium_ele->GetBinContent(nbin)));  
    //absDiff_factorized_scalefactors_Tight_ele  ->SetBinContent(nbin, TMath::Abs(absDiff_factorized_scalefactors_Tight_ele ->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Veto_ele   ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Veto_ele  ->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Loose_ele  ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Loose_ele ->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Medium_ele ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Medium_ele->GetBinContent(nbin)));  
    absDiff_unfactorized_scalefactors_Tight_ele  ->SetBinContent(nbin, TMath::Abs(absDiff_unfactorized_scalefactors_Tight_ele ->GetBinContent(nbin)));  
  }}

  // Start drawing stuff 
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;
  TLine *oneline=new TLine(10,1,100,1);
  oneline->SetLineColor(1);
  oneline->SetLineStyle(3);
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

  //TCanvas *c_absDiff_factorized_scalefactors_Veto_ele = new TCanvas("c_absDiff_factorized_scalefactors_Veto_ele","Method |difference| for Veto electron scale factors (factorized)",800,800);
  //absDiff_factorized_scalefactors_Veto_ele->Draw("TEXTE COLZ");
  //gPad->Update(); 
  //absDiff_factorized_scalefactors_Veto_ele->SetTitle("Method |difference| for Veto Electron scale factors (factorized)");
  //absDiff_factorized_scalefactors_Veto_ele->GetXaxis()->SetTitle("| #eta |");
  //absDiff_factorized_scalefactors_Veto_ele->GetXaxis()->SetTitleOffset(0.9);
  //absDiff_factorized_scalefactors_Veto_ele->GetXaxis()->SetTitleSize(0.04);
  //absDiff_factorized_scalefactors_Veto_ele->GetXaxis()->SetLabelSize(0.02);
  //absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  //absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetTitleOffset(0.9);
  //absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetTitleSize(0.04);
  //absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetLabelSize(0.02);
  //absDiff_factorized_scalefactors_Veto_ele->GetYaxis()->SetRangeUser(10,100);
  //absDiff_factorized_scalefactors_Veto_ele->SetMinimum(0);
  //absDiff_factorized_scalefactors_Veto_ele->SetMaximum(0.2);
  //absDiff_factorized_scalefactors_Veto_ele->SetMarkerSize(.8);
  //palette_axis = (TPaletteAxis*) absDiff_factorized_scalefactors_Veto_ele->GetListOfFunctions()->FindObject("palette"); 
  //palette_axis->SetLabelSize(0.02);
  //c_absDiff_factorized_scalefactors_Veto_ele->Update();
  //c_absDiff_factorized_scalefactors_Veto_ele->Print((plots_dir_method1+"absDiff_factorized_scalefactors_Veto_ele.png").c_str());
  //c_absDiff_factorized_scalefactors_Veto_ele->Print((plots_dir_method2+"absDiff_factorized_scalefactors_Veto_ele.png").c_str());

  //TCanvas *c_absDiff_factorized_scalefactors_Loose_ele = new TCanvas("c_absDiff_factorized_scalefactors_Loose_ele","Method |difference| for Loose electron scale factors (factorized)",800,800);
  //absDiff_factorized_scalefactors_Loose_ele->Draw("TEXTE COLZ");
  //gPad->Update(); 
  //absDiff_factorized_scalefactors_Loose_ele->SetTitle("Method |difference| for Loose Electron scale factors (factorized)");
  //absDiff_factorized_scalefactors_Loose_ele->GetXaxis()->SetTitle("| #eta |");
  //absDiff_factorized_scalefactors_Loose_ele->GetXaxis()->SetTitleOffset(0.9);
  //absDiff_factorized_scalefactors_Loose_ele->GetXaxis()->SetTitleSize(0.04);
  //absDiff_factorized_scalefactors_Loose_ele->GetXaxis()->SetLabelSize(0.02);
  //absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  //absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetTitleOffset(0.9);
  //absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetTitleSize(0.04);
  //absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetLabelSize(0.02);
  //absDiff_factorized_scalefactors_Loose_ele->GetYaxis()->SetRangeUser(10,100);
  //absDiff_factorized_scalefactors_Loose_ele->SetMinimum(0);
  //absDiff_factorized_scalefactors_Loose_ele->SetMaximum(0.2);
  //absDiff_factorized_scalefactors_Loose_ele->SetMarkerSize(.8);
  //palette_axis = (TPaletteAxis*) absDiff_factorized_scalefactors_Loose_ele->GetListOfFunctions()->FindObject("palette"); 
  //palette_axis->SetLabelSize(0.02);
  //c_absDiff_factorized_scalefactors_Loose_ele->Update();
  //c_absDiff_factorized_scalefactors_Loose_ele->Print((plots_dir_method1+"absDiff_factorized_scalefactors_Loose_ele.png").c_str());
  //c_absDiff_factorized_scalefactors_Loose_ele->Print((plots_dir_method2+"absDiff_factorized_scalefactors_Loose_ele.png").c_str());

  //TCanvas *c_absDiff_factorized_scalefactors_Medium_ele = new TCanvas("c_absDiff_factorized_scalefactors_Medium_ele","Method |difference| for Medium electron scale factors (factorized)",800,800);
  //absDiff_factorized_scalefactors_Medium_ele->Draw("TEXTE COLZ");
  //gPad->Update(); 
  //absDiff_factorized_scalefactors_Medium_ele->SetTitle("Method |difference| for Medium Electron scale factors (factorized)");
  //absDiff_factorized_scalefactors_Medium_ele->GetXaxis()->SetTitle("| #eta |");
  //absDiff_factorized_scalefactors_Medium_ele->GetXaxis()->SetTitleOffset(0.9);
  //absDiff_factorized_scalefactors_Medium_ele->GetXaxis()->SetTitleSize(0.04);
  //absDiff_factorized_scalefactors_Medium_ele->GetXaxis()->SetLabelSize(0.02);
  //absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  //absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetTitleOffset(0.9);
  //absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetTitleSize(0.04);
  //absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetLabelSize(0.02);
  //absDiff_factorized_scalefactors_Medium_ele->GetYaxis()->SetRangeUser(10,100);
  //absDiff_factorized_scalefactors_Medium_ele->SetMinimum(0);
  //absDiff_factorized_scalefactors_Medium_ele->SetMaximum(0.2);
  //absDiff_factorized_scalefactors_Medium_ele->SetMarkerSize(.8);
  //palette_axis = (TPaletteAxis*) absDiff_factorized_scalefactors_Medium_ele->GetListOfFunctions()->FindObject("palette"); 
  //palette_axis->SetLabelSize(0.02);
  //c_absDiff_factorized_scalefactors_Medium_ele->Update();
  //c_absDiff_factorized_scalefactors_Medium_ele->Print((plots_dir_method1+"absDiff_factorized_scalefactors_Medium_ele.png").c_str());
  //c_absDiff_factorized_scalefactors_Medium_ele->Print((plots_dir_method2+"absDiff_factorized_scalefactors_Medium_ele.png").c_str());

  //TCanvas *c_absDiff_factorized_scalefactors_Tight_ele = new TCanvas("c_absDiff_factorized_scalefactors_Tight_ele","Method |difference| for Tight electron scale factors (factorized)",800,800);
  //absDiff_factorized_scalefactors_Tight_ele->Draw("TEXTE COLZ");
  //gPad->Update(); 
  //absDiff_factorized_scalefactors_Tight_ele->SetTitle("Method |difference| for Tight Electron scale factors (factorized)");
  //absDiff_factorized_scalefactors_Tight_ele->GetXaxis()->SetTitle("| #eta |");
  //absDiff_factorized_scalefactors_Tight_ele->GetXaxis()->SetTitleOffset(0.9);
  //absDiff_factorized_scalefactors_Tight_ele->GetXaxis()->SetTitleSize(0.04);
  //absDiff_factorized_scalefactors_Tight_ele->GetXaxis()->SetLabelSize(0.02);
  //absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  //absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetTitleOffset(0.9);
  //absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetTitleSize(0.04);
  //absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetLabelSize(0.02);
  //absDiff_factorized_scalefactors_Tight_ele->GetYaxis()->SetRangeUser(10,100);
  //absDiff_factorized_scalefactors_Tight_ele->SetMinimum(0);
  //absDiff_factorized_scalefactors_Tight_ele->SetMaximum(0.2);
  //absDiff_factorized_scalefactors_Tight_ele->SetMarkerSize(.8);
  //palette_axis = (TPaletteAxis*) absDiff_factorized_scalefactors_Tight_ele->GetListOfFunctions()->FindObject("palette"); 
  //palette_axis->SetLabelSize(0.02);
  //c_absDiff_factorized_scalefactors_Tight_ele->Update();
  //c_absDiff_factorized_scalefactors_Tight_ele->Print((plots_dir_method1+"absDiff_factorized_scalefactors_Tight_ele.png").c_str());
  //c_absDiff_factorized_scalefactors_Tight_ele->Print((plots_dir_method2+"absDiff_factorized_scalefactors_Tight_ele.png").c_str());

  TCanvas *c_absDiff_unfactorized_scalefactors_Veto_ele = new TCanvas("c_absDiff_unfactorized_scalefactors_Veto_ele","Method |difference| for Veto electron scale factors",800,800);
  absDiff_unfactorized_scalefactors_Veto_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_unfactorized_scalefactors_Veto_ele->SetTitle("Method |difference| for Veto Electron scale factors");
  absDiff_unfactorized_scalefactors_Veto_ele->GetXaxis()->SetTitle("| #eta |");
  absDiff_unfactorized_scalefactors_Veto_ele->GetXaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Veto_ele->GetXaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Veto_ele->GetXaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Veto_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_unfactorized_scalefactors_Veto_ele->GetYaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Veto_ele->GetYaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Veto_ele->GetYaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Veto_ele->GetYaxis()->SetRangeUser(10,100);
  absDiff_unfactorized_scalefactors_Veto_ele->SetMinimum(0);
  absDiff_unfactorized_scalefactors_Veto_ele->SetMaximum(0.2);
  absDiff_unfactorized_scalefactors_Veto_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_unfactorized_scalefactors_Veto_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_unfactorized_scalefactors_Veto_ele->Update();
  c_absDiff_unfactorized_scalefactors_Veto_ele->Print((plots_dir_method1+"absDiff_unfactorized_scalefactors_Veto_ele.png").c_str());
  c_absDiff_unfactorized_scalefactors_Veto_ele->Print((plots_dir_method2+"absDiff_unfactorized_scalefactors_Veto_ele.png").c_str());

  TCanvas *c_absDiff_unfactorized_scalefactors_Loose_ele = new TCanvas("c_absDiff_unfactorized_scalefactors_Loose_ele","Method |difference| for Loose electron scale factors",800,800);
  absDiff_unfactorized_scalefactors_Loose_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_unfactorized_scalefactors_Loose_ele->SetTitle("Method |difference| for Loose Electron scale factors");
  absDiff_unfactorized_scalefactors_Loose_ele->GetXaxis()->SetTitle("| #eta |");
  absDiff_unfactorized_scalefactors_Loose_ele->GetXaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Loose_ele->GetXaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Loose_ele->GetXaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Loose_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_unfactorized_scalefactors_Loose_ele->GetYaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Loose_ele->GetYaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Loose_ele->GetYaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Loose_ele->GetYaxis()->SetRangeUser(10,100);
  absDiff_unfactorized_scalefactors_Loose_ele->SetMinimum(0);
  absDiff_unfactorized_scalefactors_Loose_ele->SetMaximum(0.2);
  absDiff_unfactorized_scalefactors_Loose_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_unfactorized_scalefactors_Loose_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_unfactorized_scalefactors_Loose_ele->Update();
  c_absDiff_unfactorized_scalefactors_Loose_ele->Print((plots_dir_method1+"absDiff_unfactorized_scalefactors_Loose_ele.png").c_str());
  c_absDiff_unfactorized_scalefactors_Loose_ele->Print((plots_dir_method2+"absDiff_unfactorized_scalefactors_Loose_ele.png").c_str());

  TCanvas *c_absDiff_unfactorized_scalefactors_Medium_ele = new TCanvas("c_absDiff_unfactorized_scalefactors_Medium_ele","Method |difference| for Medium electron scale factors",800,800);
  absDiff_unfactorized_scalefactors_Medium_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_unfactorized_scalefactors_Medium_ele->SetTitle("Method |difference| for Medium Electron scale factors");
  absDiff_unfactorized_scalefactors_Medium_ele->GetXaxis()->SetTitle("| #eta |");
  absDiff_unfactorized_scalefactors_Medium_ele->GetXaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Medium_ele->GetXaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Medium_ele->GetXaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Medium_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_unfactorized_scalefactors_Medium_ele->GetYaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Medium_ele->GetYaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Medium_ele->GetYaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Medium_ele->GetYaxis()->SetRangeUser(10,100);
  absDiff_unfactorized_scalefactors_Medium_ele->SetMinimum(0);
  absDiff_unfactorized_scalefactors_Medium_ele->SetMaximum(0.2);
  absDiff_unfactorized_scalefactors_Medium_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_unfactorized_scalefactors_Medium_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_unfactorized_scalefactors_Medium_ele->Update();
  c_absDiff_unfactorized_scalefactors_Medium_ele->Print((plots_dir_method1+"absDiff_unfactorized_scalefactors_Medium_ele.png").c_str());
  c_absDiff_unfactorized_scalefactors_Medium_ele->Print((plots_dir_method2+"absDiff_unfactorized_scalefactors_Medium_ele.png").c_str());

  TCanvas *c_absDiff_unfactorized_scalefactors_Tight_ele = new TCanvas("c_absDiff_unfactorized_scalefactors_Tight_ele","Method |difference| for Tight electron scale factors",800,800);
  absDiff_unfactorized_scalefactors_Tight_ele->Draw("TEXTE COLZ");
  gPad->Update(); 
  absDiff_unfactorized_scalefactors_Tight_ele->SetTitle("Method |difference| for Tight Electron scale factors");
  absDiff_unfactorized_scalefactors_Tight_ele->GetXaxis()->SetTitle("| #eta |");
  absDiff_unfactorized_scalefactors_Tight_ele->GetXaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Tight_ele->GetXaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Tight_ele->GetXaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Tight_ele->GetYaxis()->SetTitle("p_{T} [GeV]");
  absDiff_unfactorized_scalefactors_Tight_ele->GetYaxis()->SetTitleOffset(0.9);
  absDiff_unfactorized_scalefactors_Tight_ele->GetYaxis()->SetTitleSize(0.04);
  absDiff_unfactorized_scalefactors_Tight_ele->GetYaxis()->SetLabelSize(0.02);
  absDiff_unfactorized_scalefactors_Tight_ele->GetYaxis()->SetRangeUser(10,100);
  absDiff_unfactorized_scalefactors_Tight_ele->SetMinimum(0);
  absDiff_unfactorized_scalefactors_Tight_ele->SetMaximum(0.2);
  absDiff_unfactorized_scalefactors_Tight_ele->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) absDiff_unfactorized_scalefactors_Tight_ele->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_absDiff_unfactorized_scalefactors_Tight_ele->Update();
  c_absDiff_unfactorized_scalefactors_Tight_ele->Print((plots_dir_method1+"absDiff_unfactorized_scalefactors_Tight_ele.png").c_str());
  c_absDiff_unfactorized_scalefactors_Tight_ele->Print((plots_dir_method2+"absDiff_unfactorized_scalefactors_Tight_ele.png").c_str());

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
  //ele_scalefactors_file_method1->Delete("absDiff_factorized_scalefactors_Veto_ele;*");
  //ele_scalefactors_file_method1->Delete("absDiff_factorized_scalefactors_Loose_ele;*");
  //ele_scalefactors_file_method1->Delete("absDiff_factorized_scalefactors_Medium_ele;*");
  //ele_scalefactors_file_method1->Delete("absDiff_factorized_scalefactors_Tight_ele;*");
  ele_scalefactors_file_method1->Delete("absDiff_unfactorized_scalefactors_Veto_ele;*");
  ele_scalefactors_file_method1->Delete("absDiff_unfactorized_scalefactors_Loose_ele;*");
  ele_scalefactors_file_method1->Delete("absDiff_unfactorized_scalefactors_Medium_ele;*");
  ele_scalefactors_file_method1->Delete("absDiff_unfactorized_scalefactors_Tight_ele;*");
  //absDiff_factorized_scalefactors_Veto_ele->Write("absDiff_factorized_scalefactors_Veto_ele");
  //absDiff_factorized_scalefactors_Loose_ele->Write("absDiff_factorized_scalefactors_Loose_ele");
  //absDiff_factorized_scalefactors_Medium_ele->Write("absDiff_factorized_scalefactors_Medium_ele");
  //absDiff_factorized_scalefactors_Tight_ele->Write("absDiff_factorized_scalefactors_Tight_ele");
  absDiff_unfactorized_scalefactors_Veto_ele->Write("absDiff_unfactorized_scalefactors_Veto_ele");
  absDiff_unfactorized_scalefactors_Loose_ele->Write("absDiff_unfactorized_scalefactors_Loose_ele");
  absDiff_unfactorized_scalefactors_Medium_ele->Write("absDiff_unfactorized_scalefactors_Medium_ele");
  absDiff_unfactorized_scalefactors_Tight_ele->Write("absDiff_unfactorized_scalefactors_Tight_ele");
  ele_scalefactors_file_method2->cd();
  //ele_scalefactors_file_method2->Delete("absDiff_factorized_scalefactors_Veto_ele;*");
  //ele_scalefactors_file_method2->Delete("absDiff_factorized_scalefactors_Loose_ele;*");
  //ele_scalefactors_file_method2->Delete("absDiff_factorized_scalefactors_Medium_ele;*");
  //ele_scalefactors_file_method2->Delete("absDiff_factorized_scalefactors_Tight_ele;*");
  ele_scalefactors_file_method2->Delete("absDiff_unfactorized_scalefactors_Veto_ele;*");
  ele_scalefactors_file_method2->Delete("absDiff_unfactorized_scalefactors_Loose_ele;*");
  ele_scalefactors_file_method2->Delete("absDiff_unfactorized_scalefactors_Medium_ele;*");
  ele_scalefactors_file_method2->Delete("absDiff_unfactorized_scalefactors_Tight_ele;*");
  //absDiff_factorized_scalefactors_Veto_ele->Write("absDiff_factorized_scalefactors_Veto_ele");
  //absDiff_factorized_scalefactors_Loose_ele->Write("absDiff_factorized_scalefactors_Loose_ele");
  //absDiff_factorized_scalefactors_Medium_ele->Write("absDiff_factorized_scalefactors_Medium_ele");
  //absDiff_factorized_scalefactors_Tight_ele->Write("absDiff_factorized_scalefactors_Tight_ele");
  absDiff_unfactorized_scalefactors_Veto_ele->Write("absDiff_unfactorized_scalefactors_Veto_ele");
  absDiff_unfactorized_scalefactors_Loose_ele->Write("absDiff_unfactorized_scalefactors_Loose_ele");
  absDiff_unfactorized_scalefactors_Medium_ele->Write("absDiff_unfactorized_scalefactors_Medium_ele");
  absDiff_unfactorized_scalefactors_Tight_ele->Write("absDiff_unfactorized_scalefactors_Tight_ele");
  printf("* press any key to continue * \n\n");
  std::cin.ignore();

  ele_scalefactors_file_method1->Close();
  ele_scalefactors_file_method2->Close();
  mu_scalefactors_file_method1->Close();
  mu_scalefactors_file_method2->Close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plot_scale_factors(string plots_dir, string root_dir) {
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;
  gStyle->SetOptStat(0);
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);
  
  TFile *f_ele_sf_template_erfcexp = TFile::Open("~/leptonScaleFactors/root/16-02-2016/template_erfcexp/scalefactors_ele.root","READ");
  
  TH2D *ele_sf_veto        = (TH2D*) f_ele_sf_template_erfcexp ->Get("unfactorized_scalefactors_Veto_ele"); 
  TH2D *ele_sf_loose       = (TH2D*) f_ele_sf_template_erfcexp ->Get("unfactorized_scalefactors_Loose_ele"); 
  TH2D *ele_sf_medium      = (TH2D*) f_ele_sf_template_erfcexp ->Get("unfactorized_scalefactors_Medium_ele"); 
  TH2D *ele_sf_tight       = (TH2D*) f_ele_sf_template_erfcexp ->Get("unfactorized_scalefactors_Tight_ele"); 
  //TH2D *mu_sf_veto         = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Veto_mu;1");
  //TH2D *mu_sf_loose        = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Loose_mu;1");
  //TH2D *mu_sf_medium       = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Medium_mu;1");
  //TH2D *mu_sf_tight        = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Tight_mu;1");
  
  char histo_name_veto[128];
  char histo_name_loose[128];
  char histo_name_medium[128];
  char histo_name_tight[128];
    char legend_label[128];
  
  TClonesArray h_ptslice_ele_sf_veto_ ("TH1D",   n_ele_eta_bins);
  TClonesArray h_ptslice_ele_sf_loose_ ("TH1D",  n_ele_eta_bins);
  TClonesArray h_ptslice_ele_sf_medium_ ("TH1D", n_ele_eta_bins);
  TClonesArray h_ptslice_ele_sf_tight_ ("TH1D",  n_ele_eta_bins);
  TLegend *legend_ptslice_ele_sf = new TLegend(0.7, 0.6, 0.9, 0.85);
  legend_ptslice_ele_sf->SetFillColor(0);
  TClonesArray h_etaslice_ele_sf_veto_ ("TH1D",   n_ele_eta_bins);
  TClonesArray h_etaslice_ele_sf_loose_ ("TH1D",  n_ele_eta_bins);
  TClonesArray h_etaslice_ele_sf_medium_ ("TH1D", n_ele_eta_bins);
  TClonesArray h_etaslice_ele_sf_tight_ ("TH1D",  n_ele_eta_bins);
  TLegend *legend_etaslice_ele_sf = new TLegend(0.6, 0.6, 0.9, 0.85);
  legend_etaslice_ele_sf->SetFillColor(0);
  for(int i=0; i<n_ele_pt_bins; i++) {
    sprintf(histo_name_veto, "h_pt%d_ele_sf_veto", i);
    sprintf(histo_name_loose, "h_pt%d_ele_sf_loose", i);
    sprintf(histo_name_medium, "h_pt%d_ele_sf_medium", i);
    sprintf(histo_name_tight, "h_pt%d_ele_sf_tight", i);
    new(h_ptslice_ele_sf_veto_[i])   TH1D(histo_name_veto, "Electron veto S.F. as #eta in p_{T} slices", n_ele_eta_bins, ele_eta_bins);
    new(h_ptslice_ele_sf_loose_[i])  TH1D(histo_name_loose, "Electron loose S.F. as #eta in p_{T} slices", n_ele_eta_bins, ele_eta_bins);
    new(h_ptslice_ele_sf_medium_[i]) TH1D(histo_name_medium, "Electron medium S.F. as #eta in p_{T} slices", n_ele_eta_bins, ele_eta_bins);
    new(h_ptslice_ele_sf_tight_[i])  TH1D(histo_name_tight, "Electron tight S.F. as #eta in p_{T} slices", n_ele_eta_bins, ele_eta_bins);
    for(int j=0; j<n_ele_eta_bins; j++) {
    
      ((TH1D*)h_ptslice_ele_sf_veto_[i])   ->SetBinContent(j+1, ele_sf_veto->GetBinContent(ele_sf_veto->GetBin(j+1, i+1)));
      ((TH1D*)h_ptslice_ele_sf_loose_[i])  ->SetBinContent(j+1, ele_sf_loose->GetBinContent(ele_sf_loose->GetBin(j+1, i+1)));
      ((TH1D*)h_ptslice_ele_sf_medium_[i]) ->SetBinContent(j+1, ele_sf_medium->GetBinContent(ele_sf_medium->GetBin(j+1, i+1)));
      ((TH1D*)h_ptslice_ele_sf_tight_[i])  ->SetBinContent(j+1, ele_sf_tight->GetBinContent(ele_sf_tight->GetBin(j+1, i+1)));
      ((TH1D*)h_ptslice_ele_sf_veto_[i])   ->SetBinError(j+1, ele_sf_veto->GetBinError(ele_sf_veto->GetBin(j+1, i+1)));
      ((TH1D*)h_ptslice_ele_sf_loose_[i])  ->SetBinError(j+1, ele_sf_loose->GetBinError(ele_sf_loose->GetBin(j+1, i+1)));
      ((TH1D*)h_ptslice_ele_sf_medium_[i]) ->SetBinError(j+1, ele_sf_medium->GetBinError(ele_sf_medium->GetBin(j+1, i+1)));
      ((TH1D*)h_ptslice_ele_sf_tight_[i])  ->SetBinError(j+1, ele_sf_tight->GetBinError(ele_sf_tight->GetBin(j+1, i+1)));

    }

    ((TH1D*)h_ptslice_ele_sf_veto_[i])->SetMarkerColor( marker_colors[i] );
    ((TH1D*)h_ptslice_ele_sf_loose_[i])->SetMarkerColor( marker_colors[i] );
    ((TH1D*)h_ptslice_ele_sf_medium_[i])->SetMarkerColor( marker_colors[i] );
    ((TH1D*)h_ptslice_ele_sf_tight_[i])->SetMarkerColor( marker_colors[i] );
    
    ((TH1D*)h_ptslice_ele_sf_veto_[i])->SetLineColor( marker_colors[i] );
    ((TH1D*)h_ptslice_ele_sf_loose_[i])->SetLineColor( marker_colors[i] );
    ((TH1D*)h_ptslice_ele_sf_medium_[i])->SetLineColor( marker_colors[i] );
    ((TH1D*)h_ptslice_ele_sf_tight_[i])->SetLineColor( marker_colors[i] );
    
    ((TH1D*)h_ptslice_ele_sf_veto_[i])->SetLineWidth(   2 );
    ((TH1D*)h_ptslice_ele_sf_loose_[i])->SetLineWidth(  2 );
    ((TH1D*)h_ptslice_ele_sf_medium_[i])->SetLineWidth( 2 );
    ((TH1D*)h_ptslice_ele_sf_tight_[i])->SetLineWidth(  2 );
    
    ((TH1D*)h_ptslice_ele_sf_veto_[i])->SetMarkerStyle( marker_styles[i] );
    ((TH1D*)h_ptslice_ele_sf_loose_[i])->SetMarkerStyle( marker_styles[i] );
    ((TH1D*)h_ptslice_ele_sf_medium_[i])->SetMarkerStyle( marker_styles[i] );
    ((TH1D*)h_ptslice_ele_sf_tight_[i])->SetMarkerStyle( marker_styles[i] );
    sprintf(legend_label, "%i < p_{T} < %i", (int) ele_pt_bins[i], (int) ele_pt_bins[i+1]);
    legend_ptslice_ele_sf->AddEntry( ((TH1D*)h_ptslice_ele_sf_veto_[i]), legend_label, "lp");
  }
  for(int i=0; i<n_ele_eta_bins; i++) {
    sprintf(histo_name_veto, "h_eta%d_ele_sf_veto", i);
    sprintf(histo_name_loose, "h_eta%d_ele_sf_loose", i);
    sprintf(histo_name_medium, "h_eta%d_ele_sf_medium", i);
    sprintf(histo_name_tight, "h_eta%d_ele_sf_tight", i);
    new(h_etaslice_ele_sf_veto_[i])   TH1D(histo_name_veto, "Electron veto S.F. as p_{T} in #eta slices", n_ele_pt_bins, ele_pt_bins);
    new(h_etaslice_ele_sf_loose_[i])  TH1D(histo_name_loose, "Electron loose S.F. as p_{T} in #eta slices", n_ele_pt_bins, ele_pt_bins);
    new(h_etaslice_ele_sf_medium_[i]) TH1D(histo_name_medium, "Electron medium S.F. as p_{T} in #eta slices", n_ele_pt_bins, ele_pt_bins);
    new(h_etaslice_ele_sf_tight_[i])  TH1D(histo_name_tight, "Electron tight S.F. as p_{T} in #eta slices", n_ele_pt_bins, ele_pt_bins);
    for(int j=0; j<n_ele_pt_bins; j++) {
    
      ((TH1D*)h_etaslice_ele_sf_veto_[i])   ->SetBinContent(j+1, ele_sf_veto->GetBinContent(ele_sf_veto->GetBin(i+1, j+1)));
      ((TH1D*)h_etaslice_ele_sf_loose_[i])  ->SetBinContent(j+1, ele_sf_loose->GetBinContent(ele_sf_loose->GetBin(i+1, j+1)));
      ((TH1D*)h_etaslice_ele_sf_medium_[i]) ->SetBinContent(j+1, ele_sf_medium->GetBinContent(ele_sf_medium->GetBin(i+1, j+1)));
      ((TH1D*)h_etaslice_ele_sf_tight_[i])  ->SetBinContent(j+1, ele_sf_tight->GetBinContent(ele_sf_tight->GetBin(i+1, j+1)));
      ((TH1D*)h_etaslice_ele_sf_veto_[i])   ->SetBinError(j+1, ele_sf_veto->GetBinError(ele_sf_veto->GetBin(i+1, j+1)));
      ((TH1D*)h_etaslice_ele_sf_loose_[i])  ->SetBinError(j+1, ele_sf_loose->GetBinError(ele_sf_loose->GetBin(i+1, j+1)));
      ((TH1D*)h_etaslice_ele_sf_medium_[i]) ->SetBinError(j+1, ele_sf_medium->GetBinError(ele_sf_medium->GetBin(i+1, j+1)));
      ((TH1D*)h_etaslice_ele_sf_tight_[i])  ->SetBinError(j+1, ele_sf_tight->GetBinError(ele_sf_tight->GetBin(i+1, j+1)));

    }

    ((TH1D*)h_etaslice_ele_sf_veto_[i])->SetMarkerColor( marker_colors[i] );
    ((TH1D*)h_etaslice_ele_sf_loose_[i])->SetMarkerColor( marker_colors[i] );
    ((TH1D*)h_etaslice_ele_sf_medium_[i])->SetMarkerColor( marker_colors[i] );
    ((TH1D*)h_etaslice_ele_sf_tight_[i])->SetMarkerColor( marker_colors[i] );
    
    ((TH1D*)h_etaslice_ele_sf_veto_[i])->SetLineColor( marker_colors[i] );
    ((TH1D*)h_etaslice_ele_sf_loose_[i])->SetLineColor( marker_colors[i] );
    ((TH1D*)h_etaslice_ele_sf_medium_[i])->SetLineColor( marker_colors[i] );
    ((TH1D*)h_etaslice_ele_sf_tight_[i])->SetLineColor( marker_colors[i] );
    
    ((TH1D*)h_etaslice_ele_sf_veto_[i])->SetLineWidth(   2 );
    ((TH1D*)h_etaslice_ele_sf_loose_[i])->SetLineWidth(  2 );
    ((TH1D*)h_etaslice_ele_sf_medium_[i])->SetLineWidth( 2 );
    ((TH1D*)h_etaslice_ele_sf_tight_[i])->SetLineWidth(  2 );
    
    ((TH1D*)h_etaslice_ele_sf_veto_[i])->SetMarkerStyle( marker_styles[i] );
    ((TH1D*)h_etaslice_ele_sf_loose_[i])->SetMarkerStyle( marker_styles[i] );
    ((TH1D*)h_etaslice_ele_sf_medium_[i])->SetMarkerStyle( marker_styles[i] );
    ((TH1D*)h_etaslice_ele_sf_tight_[i])->SetMarkerStyle( marker_styles[i] );
    sprintf(legend_label, "%.4f < | #eta | < %.4f", ele_eta_bins[i], ele_eta_bins[i+1]);
    legend_etaslice_ele_sf->AddEntry( ((TH1D*)h_etaslice_ele_sf_veto_[i]), legend_label, "lp");
  }
  TLine *oneline = new TLine(0, 1, 2.5, 1);
  oneline->SetLineStyle(2);
  TLine *oneline2 = new TLine(10, 1, 200, 1);
  oneline2->SetLineStyle(2);
  
  TCanvas *c_ptslice_ele_sf_veto = new TCanvas("c_ptslice_ele_sf_veto", "c_ptslice_ele_sf_veto");
  c_ptslice_ele_sf_veto->SetRightMargin(0.04);
  ((TH1D*)h_ptslice_ele_sf_veto_[0])->Draw("P E1 L");
  ((TH1D*)h_ptslice_ele_sf_veto_[0])->GetXaxis()->SetTitle("| #eta |");
  ((TH1D*)h_ptslice_ele_sf_veto_[0])->SetMinimum(0.8);
  ((TH1D*)h_ptslice_ele_sf_veto_[0])->SetMaximum(1.3);
  for(int i=1; i<n_ele_pt_bins; i++) {
    ((TH1D*)h_ptslice_ele_sf_veto_[i])->Draw("P E1 L SAME");
  }
  legend_ptslice_ele_sf->Draw("SAME");
  oneline->Draw("SAME");
  c_ptslice_ele_sf_veto->Print((plots_dir+"ele_sf_veto_ptslice.png").c_str());

  TCanvas *c_ptslice_ele_sf_loose = new TCanvas("c_ptslice_ele_sf_loose", "c_ptslice_ele_sf_loose");
  c_ptslice_ele_sf_loose->SetRightMargin(0.04);
  ((TH1D*)h_ptslice_ele_sf_loose_[0])->Draw("P E1 L");
  ((TH1D*)h_ptslice_ele_sf_loose_[0])->GetXaxis()->SetTitle("| #eta |");
  ((TH1D*)h_ptslice_ele_sf_loose_[0])->SetMinimum(0.8);
  ((TH1D*)h_ptslice_ele_sf_loose_[0])->SetMaximum(1.3);
  for(int i=1; i<n_ele_pt_bins; i++) {
    ((TH1D*)h_ptslice_ele_sf_loose_[i])->Draw("P E1 L SAME");
  }
  legend_ptslice_ele_sf->Draw("SAME");
  oneline->Draw("SAME");
  c_ptslice_ele_sf_loose->Print((plots_dir+"ele_sf_loose_ptslice.png").c_str());

  TCanvas *c_ptslice_ele_sf_medium = new TCanvas("c_ptslice_ele_sf_medium", "c_ptslice_ele_sf_medium");
  c_ptslice_ele_sf_medium->SetRightMargin(0.04);
  ((TH1D*)h_ptslice_ele_sf_medium_[0])->Draw("P E1 L");
  ((TH1D*)h_ptslice_ele_sf_medium_[0])->GetXaxis()->SetTitle("| #eta |");
  ((TH1D*)h_ptslice_ele_sf_medium_[0])->SetMinimum(0.8);
  ((TH1D*)h_ptslice_ele_sf_medium_[0])->SetMaximum(1.3);
  for(int i=1; i<n_ele_pt_bins; i++) {
    ((TH1D*)h_ptslice_ele_sf_medium_[i])->Draw("P E1 L SAME");
  }
  legend_ptslice_ele_sf->Draw("SAME");
  oneline->Draw("SAME");
  c_ptslice_ele_sf_medium->Print((plots_dir+"ele_sf_medium_ptslice.png").c_str());

  TCanvas *c_ptslice_ele_sf_tight = new TCanvas("c_ptslice_ele_sf_tight", "c_ptslice_ele_sf_tight");
  c_ptslice_ele_sf_tight->SetRightMargin(0.04);
  ((TH1D*)h_ptslice_ele_sf_tight_[0])->Draw("P E1 L");
  ((TH1D*)h_ptslice_ele_sf_tight_[0])->GetXaxis()->SetTitle("| #eta |");
  ((TH1D*)h_ptslice_ele_sf_tight_[0])->SetMinimum(0.8);
  ((TH1D*)h_ptslice_ele_sf_tight_[0])->SetMaximum(1.3);
  for(int i=1; i<n_ele_pt_bins; i++) {
    ((TH1D*)h_ptslice_ele_sf_tight_[i])->Draw("P E1 L SAME");
  }
  legend_ptslice_ele_sf->Draw("SAME");
  oneline->Draw("SAME");
  c_ptslice_ele_sf_tight->Print((plots_dir+"ele_sf_tight_ptslice.png").c_str());

  TCanvas *c_etaslice_ele_sf_veto = new TCanvas("c_etaslice_ele_sf_veto", "c_etaslice_ele_sf_veto");
  c_etaslice_ele_sf_veto->SetLogx();
  c_etaslice_ele_sf_veto->SetRightMargin(0.04);
  ((TH1D*)h_etaslice_ele_sf_veto_[0])->Draw("P E1 L");
  ((TH1D*)h_etaslice_ele_sf_veto_[0])->GetXaxis()->SetTitle("p_{T} [GeV]");
  ((TH1D*)h_etaslice_ele_sf_veto_[0])->SetMinimum(0.8);
  ((TH1D*)h_etaslice_ele_sf_veto_[0])->SetMaximum(1.3);
  for(int i=1; i<n_ele_pt_bins; i++) {
    ((TH1D*)h_etaslice_ele_sf_veto_[i])->Draw("P E1 L SAME");
  }
  legend_etaslice_ele_sf->Draw("SAME");
  oneline2->Draw("SAME");
  c_etaslice_ele_sf_veto->Print((plots_dir+"ele_sf_veto_etaslice.png").c_str());

  TCanvas *c_etaslice_ele_sf_loose = new TCanvas("c_etaslice_ele_sf_loose", "c_etaslice_ele_sf_loose");
  c_etaslice_ele_sf_loose->SetLogx();
  c_etaslice_ele_sf_loose->SetRightMargin(0.04);
  ((TH1D*)h_etaslice_ele_sf_loose_[0])->Draw("P E1 L");
  ((TH1D*)h_etaslice_ele_sf_loose_[0])->GetXaxis()->SetTitle("p_{T} [GeV]");
  ((TH1D*)h_etaslice_ele_sf_loose_[0])->SetMinimum(0.8);
  ((TH1D*)h_etaslice_ele_sf_loose_[0])->SetMaximum(1.3);
  for(int i=1; i<n_ele_pt_bins; i++) {
    ((TH1D*)h_etaslice_ele_sf_loose_[i])->Draw("P E1 L SAME");
  }
  legend_etaslice_ele_sf->Draw("SAME");
  oneline2->Draw("SAME");
  c_etaslice_ele_sf_loose->Print((plots_dir+"ele_sf_loose_etaslice.png").c_str());

  TCanvas *c_etaslice_ele_sf_medium = new TCanvas("c_etaslice_ele_sf_medium", "c_etaslice_ele_sf_medium");
  c_etaslice_ele_sf_medium->SetLogx();
  c_etaslice_ele_sf_medium->SetRightMargin(0.04);
  ((TH1D*)h_etaslice_ele_sf_medium_[0])->Draw("P E1 L");
  ((TH1D*)h_etaslice_ele_sf_medium_[0])->GetXaxis()->SetTitle("p_{T} [GeV]");
  ((TH1D*)h_etaslice_ele_sf_medium_[0])->SetMinimum(0.8);
  ((TH1D*)h_etaslice_ele_sf_medium_[0])->SetMaximum(1.3);
  for(int i=1; i<n_ele_pt_bins; i++) {
    ((TH1D*)h_etaslice_ele_sf_medium_[i])->Draw("P E1 L SAME");
  }
  legend_etaslice_ele_sf->Draw("SAME");
  oneline2->Draw("SAME");
  c_etaslice_ele_sf_medium->Print((plots_dir+"ele_sf_medium_etaslice.png").c_str());

  TCanvas *c_etaslice_ele_sf_tight = new TCanvas("c_etaslice_ele_sf_tight", "c_etaslice_ele_sf_tight");
  c_etaslice_ele_sf_tight->SetLogx();
  c_etaslice_ele_sf_tight->SetRightMargin(0.04);
  ((TH1D*)h_etaslice_ele_sf_tight_[0])->Draw("P E1 L");
  ((TH1D*)h_etaslice_ele_sf_tight_[0])->GetXaxis()->SetTitle("p_{T} [GeV]");
  ((TH1D*)h_etaslice_ele_sf_tight_[0])->SetMinimum(0.8);
  ((TH1D*)h_etaslice_ele_sf_tight_[0])->SetMaximum(1.3);
  for(int i=1; i<n_ele_pt_bins; i++) {
    ((TH1D*)h_etaslice_ele_sf_tight_[i])->Draw("P E1 L SAME");
  }
  legend_etaslice_ele_sf->Draw("SAME");
  oneline2->Draw("SAME");
  c_etaslice_ele_sf_tight->Print((plots_dir+"ele_sf_tight_etaslice.png").c_str());


  // Draw unfactorized scalefactor plots
  TCanvas *c_ele_sf_veto = new TCanvas("c_ele_sf_veto","Veto electron scale factors (Data/MC)",800,800);
  ele_sf_veto->Draw("TEXTE COLZ");
  gPad->Update(); 
  ele_sf_veto->SetTitle("Veto electron scale factors (Data/MC)");
  ele_sf_veto->GetXaxis()->SetTitle("| #eta |");
  ele_sf_veto->GetXaxis()->SetTitleOffset(0.9);
  ele_sf_veto->GetXaxis()->SetTitleSize(0.04);
  ele_sf_veto->GetXaxis()->SetLabelSize(0.02);
  ele_sf_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
  ele_sf_veto->GetYaxis()->SetTitleOffset(0.9);
  ele_sf_veto->GetYaxis()->SetTitleSize(0.04);
  ele_sf_veto->GetYaxis()->SetLabelSize(0.02);
  ele_sf_veto->GetYaxis()->SetRangeUser(10,100);
  ele_sf_veto->SetMinimum(0.7);
  ele_sf_veto->SetMaximum(1.3);
  ele_sf_veto->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) ele_sf_veto->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_ele_sf_veto->Update();
  c_ele_sf_veto->Print((plots_dir+"unfactorized_scalefactors_Veto_ele.png").c_str());
  TCanvas *c_ele_sf_loose = new TCanvas("c_ele_sf_loose","Loose electron scale factors (Data/MC)",800,800);
  ele_sf_loose->Draw("TEXTE COLZ");
  gPad->Update(); 
  ele_sf_loose->SetTitle("Loose electron scale factors (Data/MC)");
  ele_sf_loose->GetXaxis()->SetTitle("| #eta |");
  ele_sf_loose->GetXaxis()->SetTitleOffset(0.9);
  ele_sf_loose->GetXaxis()->SetTitleSize(0.04);
  ele_sf_loose->GetXaxis()->SetLabelSize(0.02);
  ele_sf_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
  ele_sf_loose->GetYaxis()->SetTitleOffset(0.9);
  ele_sf_loose->GetYaxis()->SetTitleSize(0.04);
  ele_sf_loose->GetYaxis()->SetLabelSize(0.02);
  ele_sf_loose->GetYaxis()->SetRangeUser(10,100);
  ele_sf_loose->SetMinimum(0.7);
  ele_sf_loose->SetMaximum(1.3);
  ele_sf_loose->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) ele_sf_loose->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_ele_sf_loose->Update();
  c_ele_sf_loose->Print((plots_dir+"unfactorized_scalefactors_Loose_ele.png").c_str());
  TCanvas *c_ele_sf_medium = new TCanvas("c_ele_sf_medium","Medium electron scale factors (Data/MC)",800,800);
  ele_sf_medium->Draw("TEXTE COLZ");
  gPad->Update(); 
  ele_sf_medium->SetTitle("Medium electron scale factors (Data/MC)");
  ele_sf_medium->GetXaxis()->SetTitle("| #eta |");
  ele_sf_medium->GetXaxis()->SetTitleOffset(0.9);
  ele_sf_medium->GetXaxis()->SetTitleSize(0.04);
  ele_sf_medium->GetXaxis()->SetLabelSize(0.02);
  ele_sf_medium->GetYaxis()->SetTitle("p_{T} [GeV]");
  ele_sf_medium->GetYaxis()->SetTitleOffset(0.9);
  ele_sf_medium->GetYaxis()->SetTitleSize(0.04);
  ele_sf_medium->GetYaxis()->SetLabelSize(0.02);
  ele_sf_medium->GetYaxis()->SetRangeUser(10,100);
  ele_sf_medium->SetMinimum(0.7);
  ele_sf_medium->SetMaximum(1.3);
  ele_sf_medium->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) ele_sf_medium->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_ele_sf_medium->Update();
  c_ele_sf_medium->Print((plots_dir+"unfactorized_scalefactors_Medium_ele.png").c_str());
  TCanvas *c_ele_sf_tight = new TCanvas("c_ele_sf_tight","Tight electron scale factors (Data/MC)",800,800);
  ele_sf_tight->Draw("TEXTE COLZ");
  gPad->Update(); 
  ele_sf_tight->SetTitle("Tight elecron scale factors (Data/MC)");
  ele_sf_tight->GetXaxis()->SetTitle("| #eta |");
  ele_sf_tight->GetXaxis()->SetTitleOffset(0.9);
  ele_sf_tight->GetXaxis()->SetTitleSize(0.04);
  ele_sf_tight->GetXaxis()->SetLabelSize(0.02);
  ele_sf_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
  ele_sf_tight->GetYaxis()->SetTitleOffset(0.9);
  ele_sf_tight->GetYaxis()->SetTitleSize(0.04);
  ele_sf_tight->GetYaxis()->SetLabelSize(0.02);
  ele_sf_tight->GetYaxis()->SetRangeUser(10,100);
  ele_sf_tight->SetMinimum(0.7);
  ele_sf_tight->SetMaximum(1.3);
  ele_sf_tight->SetMarkerSize(.8);
  palette_axis = (TPaletteAxis*) ele_sf_tight->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_ele_sf_tight->Update();
  c_ele_sf_tight->Print((plots_dir+"unfactorized_scalefactors_Tight_ele.png").c_str());
  

}
