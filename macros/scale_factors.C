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
int marker_colors[] = {1, mit_gray, mit_red, 97, 91, 8, 60};
int marker_styles[] = {20, 21, 22, 23, 33, 34, 20};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void scale_factors(string plots_dir, string root_dir, string basename_config) {
  // This function calculates the efficiencies based on output from MIT TNP saved in subdirectories of plots_dir
  // The base names of the subdirectories are recorded in the config file
  // The efficiencies are plotted in plots_dir
  // The efficiencies and scale factors are recorded in root_dir in a 
  // rootfile whose filename is taken from basename_config

  // Pad directories with a slash at the end if it's not there
  if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
  if( root_dir[root_dir.size()-1]  != '/' )  root_dir  = root_dir + "/";

  // Open the output rootfile
  string output_rootfile_name;
  size_t ext_location = basename_config.find_last_of(".");
  if(ext_location != string::npos) output_rootfile_name = "scalefactors_"+basename_config.substr(0, ext_location)+".root";
  else output_rootfile_name = "scalefactors_"+basename_config+".root";
  TFile *output_rootfile = TFile::Open((root_dir+output_rootfile_name).c_str(), "RECREATE");
  assert(output_rootfile->IsOpen() && !output_rootfile->IsZombie());

  // Now read config file. Each line has
  // <data basename> <mc basename> <selection> <flavor>
  ifstream config_stream;
  config_stream.open(basename_config.c_str());
  assert(config_stream.is_open());
  string line, data_basename, mc_basename, output_basename, selection, flavor;
  while(getline(config_stream, line)) {
    stringstream ss(line);
    ss >> data_basename >> mc_basename >> selection >> flavor;

    if(flavor == "Electron")  output_basename = selection+"_Electron";
    else if(flavor == "Muon")   output_basename = selection+"_Muon";
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
      // Choose the max for the s.f. and eff. histogram errors
      // Analyzers who naively use this uncertainty value will just get worse sensitivity :^)
      h_sf->SetBinError(nbin, TMath::Max(
        h_sf_error_hi->GetBinContent(nbin),
        h_sf_error_lo->GetBinContent(nbin)
      ));
      h_eff_data->SetBinError(nbin, TMath::Max(
        h_error_lo_data->GetBinContent(nbin),
        h_error_hi_data->GetBinContent(nbin)
      ));
      h_eff_mc->SetBinError(nbin, TMath::Max(
        h_error_lo_mc->GetBinContent(nbin),
        h_error_hi_mc->GetBinContent(nbin)
      ));
    }} 

    // Start drawing stuff 
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.3f");
    TPaletteAxis *palette_axis;

    mitPalette();
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
    h_eff_data->GetYaxis()->SetRangeUser(10, TMath::Min(h_eff_data->GetYaxis()->GetBinUpEdge(h_eff_data->GetNbinsY()), 200.));
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
    h_eff_mc->GetYaxis()->SetRangeUser(10, TMath::Min(h_eff_mc->GetYaxis()->GetBinUpEdge(h_eff_mc->GetNbinsY()), 200.));
    h_eff_mc->SetMinimum(0);
    h_eff_mc->SetMaximum(1.);
    h_eff_mc->SetMarkerSize(.9);
    palette_axis = (TPaletteAxis*) h_eff_mc->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    canvas->Update();
    canvas->Print((plots_dir + string(h_eff_mc->GetName()) + ".png").c_str());

    mitPalette2();
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
    h_sf->GetYaxis()->SetRangeUser(10, TMath::Min(h_sf->GetYaxis()->GetBinUpEdge(h_sf->GetNbinsY()), 200.));
    h_sf->SetMinimum(.5);
    h_sf->SetMaximum(1.5);
    h_sf->SetMarkerSize(.9);
    palette_axis = (TPaletteAxis*) h_sf->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    canvas->Update();
    canvas->Print((plots_dir + string(h_sf->GetName()) + ".png").c_str());

    mitPalette();
    h_sf_error_lo->Draw("TEXTE COLZ");
    canvas->Update();
    h_sf_error_lo->GetXaxis()->SetTitle("| #eta |");
    h_sf_error_lo->GetXaxis()->SetTitleOffset(0.9);
    h_sf_error_lo->GetXaxis()->SetTitleSize(0.04);
    h_sf_error_lo->GetXaxis()->SetLabelSize(0.02);
    h_sf_error_lo->GetYaxis()->SetTitle("p_{T} [GeV]");
    h_sf_error_lo->GetYaxis()->SetTitleOffset(0.9);
    h_sf_error_lo->GetYaxis()->SetTitleSize(0.04);
    h_sf_error_lo->GetYaxis()->SetLabelSize(0.02);
    h_sf_error_lo->GetYaxis()->SetRangeUser(10, TMath::Min(h_sf_error_lo->GetYaxis()->GetBinUpEdge(h_sf_error_lo->GetNbinsY()), 200.));
    h_sf_error_lo->SetMinimum(0);
    h_sf_error_lo->SetMaximum(0.5);
    h_sf_error_lo->SetMarkerSize(.9);
    palette_axis = (TPaletteAxis*) h_sf_error_lo->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    canvas->Update();
    canvas->Print((plots_dir + string(h_sf_error_lo->GetName()) + ".png").c_str());
    h_sf_error_hi->Draw("TEXTE COLZ");
    canvas->Update();
    h_sf_error_hi->GetXaxis()->SetTitle("| #eta |");
    h_sf_error_hi->GetXaxis()->SetTitleOffset(0.9);
    h_sf_error_hi->GetXaxis()->SetTitleSize(0.04);
    h_sf_error_hi->GetXaxis()->SetLabelSize(0.02);
    h_sf_error_hi->GetYaxis()->SetTitle("p_{T} [GeV]");
    h_sf_error_hi->GetYaxis()->SetTitleOffset(0.9);
    h_sf_error_hi->GetYaxis()->SetTitleSize(0.04);
    h_sf_error_hi->GetYaxis()->SetLabelSize(0.02);
    h_sf_error_hi->GetYaxis()->SetRangeUser(10, TMath::Min(h_sf_error_hi->GetYaxis()->GetBinUpEdge(h_sf_error_hi->GetNbinsY()), 200.));
    h_sf_error_hi->SetMinimum(0);
    h_sf_error_hi->SetMaximum(0.5);
    h_sf_error_hi->SetMarkerSize(.9);
    palette_axis = (TPaletteAxis*) h_sf_error_hi->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    canvas->Update();
    canvas->Print((plots_dir + string(h_sf_error_hi->GetName()) + ".png").c_str());

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
  printf("Saved efficiencies, scale factors, and statistical errors in %s\n", (root_dir+output_rootfile_name).c_str());
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void systematics(string methods_config, string basename_config) {
  // This function calculates the systematic uncertainties by taking absolute differences of scale factor methods
  // The methods_config file format is <systematic_source> <plots_dir> <root_dir>
  // (use underscores instead of spaces in <systematic_source>
  // Output from MIT TNP are saved in subdirectories of <plots_dir>
  // Need to run scale_factors function first
  //
  // The first line of methods_config is the nominal method
  // This loops over the basename_config file and for each selection within,
  // finds the differences between the nominal method and the others,
  // then updates the root file for the nominal method with systematic uncertainties
  // and makes plots in the nominal methods plot directory <plots_dir>.


  // Graphics settings 
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);

  // First read config file. Each line has
  // <data basename> <mc basename> <selection> <flavor>
  ifstream config_stream;
  config_stream.open(basename_config.c_str());
  assert(config_stream.is_open());
  string line, data_basename, mc_basename, output_basename, selection, flavor;
  while(getline(config_stream, line)) {
    stringstream ss(line);
    ss >> data_basename >> mc_basename >> selection >> flavor;

    if(flavor == "Electron")  output_basename = selection+"_Electron";
    else if(flavor == "Muon")   output_basename = selection+"_Muon";
    else                    output_basename = selection+"_"+flavor;
    
    ifstream methods_stream;
    methods_stream.open(methods_config.c_str());
    assert(methods_stream.is_open());
    string methods_line, systematic_source, plots_dir, root_dir, nominal_plots_dir, nominal_root_dir;
    getline(methods_stream, methods_line);
    ss.str(methods_line);
    ss.clear();
    ss >> systematic_source >> nominal_plots_dir >> nominal_root_dir;
    
    // Pad directories with a slash at the end if it's not there
    if( nominal_plots_dir[plots_dir.size()-1]  != '/' ) nominal_plots_dir = nominal_plots_dir + "/";
    if( nominal_root_dir[plots_dir.size()-1]  != '/' )  nominal_root_dir  = nominal_root_dir + "/";
    
    // Open the output rootfile
    string output_rootfile_name;
    size_t ext_location = basename_config.find_last_of(".");
    if(ext_location != string::npos) output_rootfile_name = "scalefactors_"+basename_config.substr(0, ext_location)+".root";
    else output_rootfile_name = "scalefactors_"+basename_config+".root";
    TFile *output_rootfile = TFile::Open((nominal_root_dir+output_rootfile_name).c_str(), "UPDATE");
    assert(output_rootfile->IsOpen() && !output_rootfile->IsZombie());

    TH2D *h_nominal_sf = (TH2D*) output_rootfile->Get(("scalefactors_"+output_basename).c_str());
    TH2D *h_syst_combined = (TH2D*) h_nominal_sf->Clone();
    h_syst_combined->Reset();
    h_syst_combined->SetName(("scalefactors_"+output_basename+"_syst_error_combined").c_str());
    h_syst_combined->SetTitle(("Combined systematics for "+selection+" "+flavor+" selection").c_str());
    
    // Now loop over the different methods and compute absolute differences in scale factors
    while(getline(methods_stream, methods_line)) {
      ss.str(methods_line);
      ss.clear();
      ss >> systematic_source >> plots_dir >> root_dir;
      if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
      if( root_dir[plots_dir.size()-1]  != '/' )  root_dir  = root_dir + "/";
      //printf("root dir is %s\n",  root_dir.c_str());
      // replace underscores in the method name for the plots
      string nice_method_name = systematic_source;
      size_t pos = 0;
      while((pos = nice_method_name.find("_", pos)) != std::string::npos){
         nice_method_name.replace(pos, 1, " ");
         pos += 1;
      } 
      TH2D *h_syst_method = (TH2D*) h_syst_combined->Clone(); 
      h_syst_method->Scale(0);
      h_syst_method->SetName(("scalefactors_"+output_basename+"_syst_error_"+systematic_source).c_str());
      h_syst_method->SetTitle(("Systematics for "+selection+" "+flavor+" selection ("+nice_method_name+")").c_str());

      // open the file with the alternate method scale factors
      TFile *alternate_rootfile = TFile::Open((root_dir+output_rootfile_name).c_str(), "READ");
      assert(alternate_rootfile->IsOpen() && !alternate_rootfile->IsZombie());

      TH2D *h_alternate_sf = (TH2D*) alternate_rootfile->Get(("scalefactors_"+output_basename).c_str());
      for(int i = 1; i <= h_nominal_sf->GetNbinsX(); i++) { for(int j = 1; j <= h_nominal_sf->GetNbinsY(); j++) {
        unsigned int nbin = h_nominal_sf->GetBin(i,j);
        h_syst_method->SetBinContent(nbin, TMath::Abs(h_nominal_sf->GetBinContent(nbin) - h_alternate_sf->GetBinContent(nbin)));
        h_syst_combined->SetBinContent(nbin, h_syst_combined->GetBinContent(nbin) + pow( h_syst_method->GetBinContent(nbin), 2));
      }}

      // Draw 2D histogram of systematics for this method 
      TCanvas *canvas = new TCanvas("canvas", "canvas", 800,600);
      h_syst_method->SetMinimum(0);
      h_syst_method->SetMaximum(.2);
      h_syst_method->Draw("TEXTE COLZ");
      canvas->Update();
      h_syst_method->GetXaxis()->SetTitle("| #eta |");
      h_syst_method->GetXaxis()->SetTitleOffset(0.9);
      h_syst_method->GetXaxis()->SetTitleSize(0.04);
      h_syst_method->GetXaxis()->SetLabelSize(0.02);
      h_syst_method->GetYaxis()->SetTitle("p_{T} [GeV]");
      h_syst_method->GetYaxis()->SetTitleOffset(0.9);
      h_syst_method->GetYaxis()->SetTitleSize(0.04);
      h_syst_method->GetYaxis()->SetLabelSize(0.02);
      h_syst_method->GetYaxis()->SetRangeUser(         10., TMath::Min(h_syst_method->GetYaxis()->GetBinUpEdge( h_syst_method->GetNbinsY() ) , 200.)       );
      h_syst_method->SetMarkerSize(.9);
      palette_axis = (TPaletteAxis*) h_syst_method->GetListOfFunctions()->FindObject("palette"); 
      palette_axis->SetLabelSize(0.02);
      canvas->Update();
      canvas->Print((nominal_plots_dir + string(h_syst_method->GetName()) + ".png").c_str());

      output_rootfile->cd();
      h_syst_method->Write();
      alternate_rootfile->Close();
      delete alternate_rootfile;
      delete h_syst_method;
      delete canvas;
    }
    // get the hi/low stat. errors
    TH2D *h_nominal_sf_stat_error_lo = (TH2D*) output_rootfile->Get(("scalefactors_"+output_basename+"_stat_error_lo").c_str());
    TH2D *h_nominal_sf_stat_error_hi = (TH2D*) output_rootfile->Get(("scalefactors_"+output_basename+"_stat_error_hi").c_str());

    TH2D *h_sf  = (TH2D*) h_nominal_sf->Clone();
    for(int i = 1; i <= h_nominal_sf->GetNbinsX(); i++) { for(int j = 1; j <= h_nominal_sf->GetNbinsY(); j++) {
      unsigned int nbin = h_nominal_sf->GetBin(i,j);
      double stat_error_lo = h_nominal_sf_stat_error_lo->GetBinContent(nbin);
      double stat_error_hi = h_nominal_sf_stat_error_hi->GetBinContent(nbin);
      h_syst_combined->SetBinContent(nbin, sqrt(h_syst_combined->GetBinContent(nbin)));
      printf("setting bin error for bin %d = %f\n", nbin, sqrt(
        pow(h_syst_combined->GetBinContent(nbin),2) +
        pow(TMath::Max(stat_error_lo, stat_error_hi), 2)
      ));
      h_sf->SetBinError(nbin, sqrt(
        pow(h_syst_combined->GetBinContent(nbin),2) +
        pow(TMath::Max(stat_error_lo, stat_error_hi), 2)
      ));
    }}
    
    // Redraw 2D SF histograms with systematics and statistical error    
    mitPalette2();
    TCanvas *c_sf = new TCanvas("c_sf", "c_sf", 800,600);
    h_sf->SetMinimum(.8);
    h_sf->SetMaximum(1.2);
    h_sf->Draw("TEXTE COLZ");
    c_sf->Update();
    h_sf->GetXaxis()->SetTitle("| #eta |");
    h_sf->GetXaxis()->SetTitleOffset(0.9);
    h_sf->GetXaxis()->SetTitleSize(0.04);
    h_sf->GetXaxis()->SetLabelSize(0.02);
    h_sf->GetYaxis()->SetTitle("p_{T} [GeV]");
    h_sf->GetYaxis()->SetTitleOffset(0.9);
    h_sf->GetYaxis()->SetTitleSize(0.04);
    h_sf->GetYaxis()->SetLabelSize(0.02);
    h_sf->GetYaxis()->SetRangeUser(
      10., TMath::Min(h_sf->GetYaxis()->GetBinUpEdge( h_sf->GetNbinsY() ) , 200.)
    );
    h_sf->SetMarkerSize(.9);
    palette_axis = (TPaletteAxis*) h_sf->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_sf->Update();
    c_sf->Print((plots_dir + string(h_nominal_sf->GetName()) + ".png").c_str());

    mitPalette();
    // Draw 2D histogram of all systematics combined in quadrature    
    TCanvas *c_syst_combined = new TCanvas("c_syst_combined", "c_syst_combined", 800,600);
    h_syst_combined->SetMinimum(0);
    h_syst_combined->SetMaximum(.2);
    h_syst_combined->Draw("TEXTE COLZ");
    c_syst_combined->Update();
    h_syst_combined->GetXaxis()->SetTitle("| #eta |");
    h_syst_combined->GetXaxis()->SetTitleOffset(0.9);
    h_syst_combined->GetXaxis()->SetTitleSize(0.04);
    h_syst_combined->GetXaxis()->SetLabelSize(0.02);
    h_syst_combined->GetYaxis()->SetTitle("p_{T} [GeV]");
    h_syst_combined->GetYaxis()->SetTitleOffset(0.9);
    h_syst_combined->GetYaxis()->SetTitleSize(0.04);
    h_syst_combined->GetYaxis()->SetLabelSize(0.02);
    h_syst_combined->GetYaxis()->SetRangeUser(
      10., TMath::Min(h_syst_combined->GetYaxis()->GetBinUpEdge( h_syst_combined->GetNbinsY() ) , 200.)
    );
    h_syst_combined->SetMarkerSize(.9);
    palette_axis = (TPaletteAxis*) h_syst_combined->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_combined->Update();
    c_syst_combined->Print((nominal_plots_dir + string(h_syst_combined->GetName()) + ".png").c_str());
    
    // Now draw 1D histograms of scale factors in eta and pT slices with full errors

    int maxslices=7; 
     
    // do pt slices
    TObjArray *ptslices = new TObjArray( h_nominal_sf->GetNbinsY() );
    ptslices->SetOwner(kTRUE);
    for(int j = 1; j <= TMath::Min(maxslices, h_nominal_sf->GetNbinsY()); j++) { // loop over each pT slice
      TGraphAsymmErrors *ptslice = new TGraphAsymmErrors( h_nominal_sf->GetNbinsX() );
      for(int i = 1; i <= TMath::Min(maxslices, h_nominal_sf->GetNbinsX()); i++) { // fill the eta points
        double bincenter = h_nominal_sf->GetXaxis()->GetBinCenter(i);
        double binwidth = h_nominal_sf->GetXaxis()->GetBinUpEdge(i) - h_nominal_sf->GetXaxis()->GetBinLowEdge(i);
        double nbin = h_nominal_sf->GetBin(i,j);
        double stat_error_lo = h_nominal_sf_stat_error_lo->GetBinContent(nbin);
        double stat_error_hi = h_nominal_sf_stat_error_hi->GetBinContent(nbin);
        double syst_error    = h_syst_combined->GetBinContent(nbin);
        ptslice->SetPoint(i-1, bincenter, h_nominal_sf->GetBinContent(nbin));
        ptslice->SetPointError(i-1,
          binwidth/2.,
          binwidth/2.,
          sqrt(pow(stat_error_lo, 2) + pow(syst_error, 2)),
          sqrt(pow(stat_error_hi, 2) + pow(syst_error, 2))
        );
      }
      char name[128];
      sprintf(name, "ptslice%d", j); 
      ptslice->SetName(name);
      ptslice->SetTitle((selection+" "+flavor+" scalefactors ( p_{T} slices )").c_str() );
      ptslice->SetMarkerColor(marker_colors[j-1]);
      ptslice->SetMarkerStyle(marker_styles[j-1]);
      ptslice->SetLineColor(marker_colors[j-1]);
      ptslice->SetLineWidth(2);
      ptslice->SetMinimum(0.8);
      ptslice->SetMaximum(1.2);
      ptslices->Add(ptslice);
    }
    
    TCanvas *canvas_ptslices = new TCanvas("canvas_ptslices", "canvas_ptslices", 800,600);
    canvas_ptslices->SetRightMargin(0.05);
    canvas_ptslices->SetLeftMargin(0.07);
    TLegend *legend_ptslices = new TLegend(.7,.7,.90,.88);
    legend_ptslices->SetFillColor(0);
    legend_ptslices->SetMargin(.5);
    for(int j = 1; j <= TMath::Min(maxslices, h_nominal_sf->GetNbinsY()); j++) { // loop over each pt slice
      if(j==1)   ( (TGraphAsymmErrors*) ptslices->At(j-1)) ->Draw("AP");
      else       ( (TGraphAsymmErrors*) ptslices->At(j-1)) ->Draw("P SAME");
      char legend_label[128];
      sprintf(legend_label, "%i < p_{T} < %i   ", (int) h_nominal_sf->GetYaxis()->GetBinLowEdge(j), (int) h_nominal_sf->GetYaxis()->GetBinUpEdge(j)); 
      legend_ptslices->AddEntry(( (TGraphAsymmErrors*) ptslices->At(j-1)), legend_label, "lp");
    }
    ( (TGraphAsymmErrors*) ptslices->At(0))->GetXaxis()->SetRangeUser(
      h_nominal_sf->GetXaxis()->GetBinLowEdge(1),
      h_nominal_sf->GetXaxis()->GetBinUpEdge( h_nominal_sf->GetNbinsX() )
    );
    ( (TGraphAsymmErrors*) ptslices->At(0))->GetXaxis()->SetTitle("|#eta|");
    ( (TGraphAsymmErrors*) ptslices->At(0))->GetXaxis()->SetTitleOffset(1.3);
    TLine *oneline_ptslices = new TLine(
      h_nominal_sf->GetXaxis()->GetBinLowEdge(1),
      1,
      h_nominal_sf->GetXaxis()->GetBinUpEdge( h_nominal_sf->GetNbinsX() ),
      1
    );
    oneline_ptslices->SetLineStyle(2);
    oneline_ptslices->Draw("SAME");
    legend_ptslices->Draw("SAME");
    canvas_ptslices->Print((nominal_plots_dir + "scalefactors_"+output_basename+"_ptslices.png").c_str());
  
    // Do eta slices
    TObjArray *etaslices = new TObjArray( h_nominal_sf->GetNbinsX() );
    etaslices->SetOwner(kTRUE);
    for(int i = 1; i <= TMath::Min(maxslices, h_nominal_sf->GetNbinsX()); i++) { // loop over each eta slice
      TGraphAsymmErrors *etaslice = new TGraphAsymmErrors( h_nominal_sf->GetNbinsY() );
      for(int j = 1; j <= TMath::Min(maxslices, h_nominal_sf->GetNbinsY()); j++) { // fill the pT points
        double bincenter = h_nominal_sf->GetYaxis()->GetBinCenter(j);
        double binwidth = h_nominal_sf->GetYaxis()->GetBinUpEdge(j) - h_nominal_sf->GetYaxis()->GetBinLowEdge(j);
        double nbin = h_nominal_sf->GetBin(i,j);
        double stat_error_lo = h_nominal_sf_stat_error_lo->GetBinContent(nbin);
        double stat_error_hi = h_nominal_sf_stat_error_hi->GetBinContent(nbin);
        double syst_error    = h_syst_combined->GetBinContent(nbin);
        etaslice->SetPoint(j-1, bincenter, h_nominal_sf->GetBinContent(nbin));
        etaslice->SetPointError(j-1,
          binwidth/2.,
          binwidth/2.,
          sqrt(pow(stat_error_lo, 2) + pow(syst_error, 2)),
          sqrt(pow(stat_error_hi, 2) + pow(syst_error, 2))
        );
      }
      char name[128];
      sprintf(name, "etaslice%d", i); 
      etaslice->SetName(name);
      etaslice->SetTitle((selection+" "+flavor+" scalefactors ( |#eta| slices )").c_str() );
      etaslice->SetMarkerColor(marker_colors[i-1]);
      etaslice->SetMarkerStyle(marker_styles[i-1]);
      etaslice->SetLineColor(marker_colors[i-1]);
      etaslice->SetLineWidth(2);
      etaslice->SetMinimum(0.8);
      etaslice->SetMaximum(1.2);
      etaslices->Add(etaslice);
    }
    TCanvas *canvas_etaslices = new TCanvas("canvas_etaslices", "canvas_etaslices", 800,600);
    canvas_etaslices->SetLogx();
    canvas_etaslices->SetRightMargin(0.05);
    canvas_etaslices->SetLeftMargin(0.07);
    TLegend *legend_etaslices = new TLegend(.7,.6,.9,.85);
    legend_etaslices->SetFillColor(0);
    for(int i = 1; i <= TMath::Min(maxslices, h_nominal_sf->GetNbinsX()); i++) { // loop over each eta slice
      if(i==1)   ( (TGraphAsymmErrors*) etaslices->At(i-1)) ->Draw("AP");
      else       ( (TGraphAsymmErrors*) etaslices->At(i-1)) ->Draw("P SAME");
      char legend_label[128];
      sprintf(legend_label, "%.4f < |#eta| < %.4f", h_nominal_sf->GetXaxis()->GetBinLowEdge(i), h_nominal_sf->GetXaxis()->GetBinUpEdge(i)); 
      legend_etaslices->AddEntry(( (TGraphAsymmErrors*) etaslices->At(i-1)), legend_label, "lp");
    }
    ( (TGraphAsymmErrors*) etaslices->At(0))->GetXaxis()->SetRangeUser(
      h_nominal_sf->GetYaxis()->GetBinLowEdge(1),
      h_nominal_sf->GetYaxis()->GetBinUpEdge( h_nominal_sf->GetNbinsY() )
    );
    ( (TGraphAsymmErrors*) etaslices->At(0))->GetXaxis()->SetTitle("p_{T} [GeV]");
    ( (TGraphAsymmErrors*) etaslices->At(0))->GetXaxis()->SetTitleOffset(1.3);
    ( (TGraphAsymmErrors*) etaslices->At(0))->GetXaxis()->SetMoreLogLabels();
    TLine *oneline_etaslices = new TLine(
      h_nominal_sf->GetYaxis()->GetBinLowEdge(1),
      1,
      h_nominal_sf->GetYaxis()->GetBinUpEdge( h_nominal_sf->GetNbinsY() ),
      1
    );
    oneline_etaslices->SetLineStyle(2);
    oneline_etaslices->Draw("SAME");
    legend_etaslices->Draw("SAME");
    canvas_etaslices->Print((nominal_plots_dir + "scalefactors_"+output_basename+"_etaslices.png").c_str());
    // Write to file
    output_rootfile->cd();
    h_syst_combined->Write();
    output_rootfile->Close();
    printf("Saved systematic errors for %s %s selection in %s\n", selection.c_str(), flavor.c_str(), (nominal_root_dir+output_rootfile_name).c_str());
    delete output_rootfile;

    //delete palette_axis;
    //delete h_syst_combined;
    //delete h_nominal_sf;
    delete c_sf;
    delete c_syst_combined;
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void wholething() {
  scale_factors("plots/2016-02-26/2016-02-26_74x_BWCBPlusVoigt_erfcexp/",    "root/2016-02-26/2016-02-26_74x_BWCBPlusVoigt_erfcexp/",   "ele_74x.cfg");
  scale_factors("plots/2016-02-26/2016-02-26_74x_template_erfcexp/",         "root/2016-02-26/2016-02-26_74x_template_erfcexp",         "ele_74x.cfg");
  scale_factors("plots/2016-02-26/2016-02-26_74x_template_exp/",             "root/2016-02-26/2016-02-26_74x_template_exp",             "ele_74x.cfg");
  scale_factors("plots/2016-02-26/2016-02-26_74x_template_erfcexp_diffTag/", "root/2016-02-26/2016-02-26_74x_template_erfcexp_diffTag", "ele_74x.cfg");
  scale_factors("plots/2016-02-26/2016-02-26_74x_template_erfcexp_LO/",      "root/2016-02-26/2016-02-26_74x_template_erfcexp_LO",      "ele_74x.cfg");
  scale_factors("plots/2016-02-26/2016-02-26_74x_BWCBPlusVoigt_erfcexp/",    "root/2016-02-26/2016-02-26_74x_BWCBPlusVoigt_erfcexp/",   "mu_74x.cfg");
  scale_factors("plots/2016-02-26/2016-02-26_74x_template_erfcexp/",         "root/2016-02-26/2016-02-26_74x_template_erfcexp",         "mu_74x.cfg");
  scale_factors("plots/2016-02-26/2016-02-26_74x_template_exp/",             "root/2016-02-26/2016-02-26_74x_template_exp",             "mu_74x.cfg");
  scale_factors("plots/2016-02-26/2016-02-26_74x_template_erfcexp_LO/",      "root/2016-02-26/2016-02-26_74x_template_erfcexp_LO",      "mu_74x.cfg");
  scale_factors("plots/2016-02-26/2016-02-26_74x_template_erfcexp_diffTag/", "root/2016-02-26/2016-02-26_74x_template_erfcexp_diffTag", "mu_74x.cfg");
  systematics("methods_74x.cfg", "ele_74x.cfg");
  systematics("methods_74x.cfg", "mu_74x.cfg");
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
/////////////////////////////////////////////////////////////////////////////
void plot_sf_1d(string data_dir, string mc_dir, bool do_pt=true, bool do_eta=true, bool do_phi=true, bool do_npv=true, bool do_etapt=true) {
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);

  if( data_dir[data_dir.size()-1]  != '/' ) data_dir = data_dir + "/";
  if( mc_dir[mc_dir.size()-1]  != '/' )  mc_dir  = mc_dir + "/";
  TFile *f_data=TFile::Open((data_dir+"eff.root").c_str());
  TFile *f_mc  =TFile::Open((mc_dir+"eff.root").c_str());
  TGraphAsymmErrors *g_data_pt, 
                    *g_data_eta,
                    *g_data_phi,
                    *g_data_npv,
                    *g_mc_pt,  
                    *g_mc_eta, 
                    *g_mc_phi, 
                    *g_mc_npv, 
                    *g_ratio_pt, 
                    *g_ratio_eta,
                    *g_ratio_phi,
                    *g_ratio_npv;
  bool logx;
  if(do_pt) {
    g_data_pt  = (TGraphAsymmErrors*) f_data->Get("grEffPt");
    g_mc_pt  = (TGraphAsymmErrors*) f_mc->Get("grEffPt");
    g_ratio_pt  = (TGraphAsymmErrors*) g_data_pt ->Clone("g_ratio_pt");
    g_data_pt ->SetTitle("");
    g_mc_pt   ->SetTitle("");
    g_ratio_pt ->SetTitle(""); 
    unsigned int np_pt  = g_data_pt ->GetN();
    for(unsigned int ip = 0; ip < np_pt; ip++) {
      Double_t x_data,y_data,x_mc,y_mc;
      g_data_pt->GetPoint(ip,x_data,y_data);
      g_mc_pt->GetPoint(ip,x_mc,y_mc);
      g_ratio_pt->SetPoint(ip, x_data, y_data/y_mc);
      g_ratio_pt->SetPointError(ip,
        g_data_pt->GetErrorXlow(ip),
        g_data_pt->GetErrorXhigh(ip),
        sqrt(pow(g_data_pt->GetErrorYlow(ip)/y_mc,2) + pow(g_mc_pt->GetErrorYhigh(ip)*y_data/(y_mc*y_mc),2)),
        sqrt(pow(g_data_pt->GetErrorYhigh(ip)/y_mc,2) + pow(g_mc_pt->GetErrorYlow(ip)*y_data/(y_mc*y_mc),2))
      );
    }
    logx=true;
    TCanvas *canvas_pt = new TCanvas("canvas_pt", "canvas_pt", 640,480);
    canvas_pt->SetMargin(0,0,0,0);
    TPad *pad1_pt = new TPad("pad1_pt", "pad1_pt", 0, .3, 1, 1);
    pad1_pt->SetGrid(0,1);
    pad1_pt->SetMargin(0.1,0.04,0,.1);
    pad1_pt->Draw();
    pad1_pt->cd();
    if(logx) pad1_pt->SetLogx();
    g_mc_pt->GetYaxis()->SetTitleSize(15);
    g_mc_pt->GetYaxis()->SetTitleFont(43);
    g_mc_pt->GetYaxis()->SetTitleOffset(1.55);
    g_mc_pt->SetMarkerColor(mit_red);
    g_mc_pt->SetLineColor(mit_red);
    g_mc_pt->GetYaxis()->SetTitle("#varepsilon");
    g_mc_pt->Draw("ap");
    g_mc_pt->SetMaximum(1.1);
    g_mc_pt->SetMinimum(0.15);
    g_data_pt->Draw("p same");
    TPad *pad2_pt = new TPad("pad2_pt", "pad2_pt", 0, 0.05, 1, 0.3);
    canvas_pt->cd();
    pad2_pt->SetMargin(0.1,0.04,0.3,0);
    pad2_pt->Draw();
    pad2_pt->cd();
    if(logx) pad2_pt->SetLogx();
    g_ratio_pt->SetMaximum(1.2);
    g_ratio_pt->SetMinimum(0.8);
    if(logx) g_ratio_pt->GetXaxis()->SetMoreLogLabels();
    g_ratio_pt->GetXaxis()->SetTitle("p_{T} [GeV]");
    g_ratio_pt->GetYaxis()->SetTitle("Data/MC");
    g_ratio_pt->GetYaxis()->SetNdivisions(5);
    g_ratio_pt->GetYaxis()->SetTitleSize(15);
    g_ratio_pt->GetYaxis()->SetTitleFont(43);
    g_ratio_pt->GetYaxis()->SetTitleOffset(1.55);
    g_ratio_pt->GetYaxis()->SetLabelFont(43); 
    g_ratio_pt->GetYaxis()->SetLabelSize(15);
    g_ratio_pt->GetXaxis()->SetTitleSize(15);
    g_ratio_pt->GetXaxis()->SetTitleFont(43);
    g_ratio_pt->GetXaxis()->SetTitleOffset(5.);
    g_ratio_pt->GetXaxis()->SetLabelFont(43);
    g_ratio_pt->GetXaxis()->SetLabelSize(15);
    TLine *oneline_pt = new TLine(g_ratio_pt->GetXaxis()->GetXmin(),1,g_ratio_pt->GetXaxis()->GetXmax(),1);
    oneline_pt->SetLineColor(mit_gray);
    oneline_pt->SetLineWidth(1);
    oneline_pt->SetLineStyle(3);
    g_ratio_pt->Draw("ap");
    oneline_pt->Draw("SAME");
    g_ratio_pt->Draw("p same");
    canvas_pt->Print((data_dir+"eff_ratio_pt.png").c_str());
  }
  if(do_eta) {
    TGraphAsymmErrors *g_data_eta = (TGraphAsymmErrors*) f_data->Get("grEffEta");
    TGraphAsymmErrors *g_mc_eta = (TGraphAsymmErrors*) f_mc->Get("grEffEta");
    TGraphAsymmErrors *g_ratio_eta = (TGraphAsymmErrors*) g_data_eta->Clone("g_ratio_eta");
    g_data_eta->SetTitle("");
    g_mc_eta  ->SetTitle("");
    g_ratio_eta->SetTitle(""); 
    unsigned int np_eta = g_data_eta->GetN();
    for(unsigned int ip = 0; ip < np_eta; ip++) {
      Double_t x_data,y_data,x_mc,y_mc;
      g_data_eta->GetPoint(ip,x_data,y_data);
      g_mc_eta->GetPoint(ip,x_mc,y_mc);
      g_ratio_eta->SetPoint(ip, x_data, y_data/y_mc);
      g_ratio_eta->SetPointError(ip,
        g_data_eta->GetErrorXlow(ip),
        g_data_eta->GetErrorXhigh(ip),
        sqrt(pow(g_data_eta->GetErrorYlow(ip)/y_mc,2) + pow(g_mc_eta->GetErrorYhigh(ip)*y_data/(y_mc*y_mc),2)),
        sqrt(pow(g_data_eta->GetErrorYhigh(ip)/y_mc,2) + pow(g_mc_eta->GetErrorYlow(ip)*y_data/(y_mc*y_mc),2))
      );
    }
    logx=false;
    TCanvas *canvas_eta = new TCanvas("canvas_eta", "canvas_eta", 640,480);
    canvas_eta->SetMargin(0,0,0,0);
    TPad *pad1_eta = new TPad("pad1_eta", "pad1_eta", 0, .3, 1, 1);
    pad1_eta->SetGrid(0,1);
    pad1_eta->SetMargin(0.1,0.04,0,.1);
    pad1_eta->Draw();
    pad1_eta->cd();
    if(logx) pad1_eta->SetLogx();
    g_mc_eta->GetYaxis()->SetTitleSize(15);
    g_mc_eta->GetYaxis()->SetTitleFont(43);
    g_mc_eta->GetYaxis()->SetTitleOffset(1.55);
    g_mc_eta->SetMarkerColor(mit_red);
    g_mc_eta->SetLineColor(mit_red);
    g_data_eta->GetYaxis()->SetTitle("#varepsilon");
    g_data_eta->Draw("ap");
    g_data_eta->SetMaximum(1.1);
    g_data_eta->SetMinimum(0.15);
    g_mc_eta->Draw("p same");
    TPad *pad2_eta = new TPad("pad2_eta", "pad2_eta", 0, 0.05, 1, 0.3);
    canvas_eta->cd();
    pad2_eta->SetMargin(0.1,0.04,0.3,0);
    pad2_eta->Draw();
    pad2_eta->cd();
    if(logx) pad2_eta->SetLogx();
    g_ratio_eta->SetMaximum(1.2);
    g_ratio_eta->SetMinimum(0.8);
    if(logx) g_ratio_eta->GetXaxis()->SetMoreLogLabels();
    g_ratio_eta->GetXaxis()->SetTitle("#eta");
    g_ratio_eta->GetYaxis()->SetTitle("Data/MC");
    g_ratio_eta->GetYaxis()->SetNdivisions(5);
    g_ratio_eta->GetYaxis()->SetTitleSize(15);
    g_ratio_eta->GetYaxis()->SetTitleFont(43);
    g_ratio_eta->GetYaxis()->SetTitleOffset(1.55);
    g_ratio_eta->GetYaxis()->SetLabelFont(43); 
    g_ratio_eta->GetYaxis()->SetLabelSize(15);
    g_ratio_eta->GetXaxis()->SetTitleSize(15);
    g_ratio_eta->GetXaxis()->SetTitleFont(43);
    g_ratio_eta->GetXaxis()->SetTitleOffset(5.);
    g_ratio_eta->GetXaxis()->SetLabelFont(43);
    g_ratio_eta->GetXaxis()->SetLabelSize(15);
    TLine *oneline_eta = new TLine(g_ratio_eta->GetXaxis()->GetXmin(),1,g_ratio_eta->GetXaxis()->GetXmax(),1);
    oneline_eta->SetLineColor(mit_gray);
    oneline_eta->SetLineWidth(1);
    oneline_eta->SetLineStyle(3);
    g_ratio_eta->Draw("ap");
    oneline_eta->Draw("SAME");
    g_ratio_eta->Draw("p same");
    canvas_eta->Print((data_dir+"eff_ratio_eta.png").c_str());
  }
  if(do_phi) {
    TGraphAsymmErrors *g_data_phi = (TGraphAsymmErrors*) f_data->Get("grEffPhi");
    TGraphAsymmErrors *g_mc_phi = (TGraphAsymmErrors*) f_mc->Get("grEffPhi");
    TGraphAsymmErrors *g_ratio_phi = (TGraphAsymmErrors*) g_data_phi->Clone("g_ratio_phi");
    g_data_phi->SetTitle("");
    g_mc_phi  ->SetTitle("");
    g_ratio_phi->SetTitle(""); 
    unsigned int np_phi = g_data_phi->GetN();
    for(unsigned int ip = 0; ip < np_phi; ip++) {
      Double_t x_data,y_data,x_mc,y_mc;
      g_data_phi->GetPoint(ip,x_data,y_data);
      g_mc_phi->GetPoint(ip,x_mc,y_mc);
      g_ratio_phi->SetPoint(ip, x_data, y_data/y_mc);
      g_ratio_phi->SetPointError(ip,
        g_data_phi->GetErrorXlow(ip),
        g_data_phi->GetErrorXhigh(ip),
        sqrt(pow(g_data_phi->GetErrorYlow(ip)/y_mc,2) + pow(g_mc_phi->GetErrorYhigh(ip)*y_data/(y_mc*y_mc),2)),
        sqrt(pow(g_data_phi->GetErrorYhigh(ip)/y_mc,2) + pow(g_mc_phi->GetErrorYlow(ip)*y_data/(y_mc*y_mc),2))
      );
    }
    logx=false;
    TCanvas *canvas_phi = new TCanvas("canvas_phi", "canvas_phi", 640,480);
    canvas_phi->SetMargin(0,0,0,0);
    TPad *pad1_phi = new TPad("pad1_phi", "pad1_phi", 0, .3, 1, 1);
    pad1_phi->SetGrid(0,1);
    pad1_phi->SetMargin(0.1,0.04,0,.1);
    pad1_phi->Draw();
    pad1_phi->cd();
    if(logx) pad1_phi->SetLogx();
    g_mc_phi->GetYaxis()->SetTitleSize(15);
    g_mc_phi->GetYaxis()->SetTitleFont(43);
    g_mc_phi->GetYaxis()->SetTitleOffset(1.55);
    g_mc_phi->SetMarkerColor(mit_red);
    g_mc_phi->SetLineColor(mit_red);
    g_data_phi->GetYaxis()->SetTitle("#varepsilon");
    g_data_phi->Draw("ap");
    g_data_phi->SetMaximum(1.1);
    g_data_phi->SetMinimum(0.15);
    g_mc_phi->Draw("p same");
    TPad *pad2_phi = new TPad("pad2_phi", "pad2_phi", 0, 0.05, 1, 0.3);
    canvas_phi->cd();
    pad2_phi->SetMargin(0.1,0.04,0.3,0);
    pad2_phi->Draw();
    pad2_phi->cd();
    if(logx) pad2_phi->SetLogx();
    g_ratio_phi->SetMaximum(1.2);
    g_ratio_phi->SetMinimum(0.8);
    if(logx) g_ratio_phi->GetXaxis()->SetMoreLogLabels();
    g_ratio_phi->GetXaxis()->SetTitle("#phi");
    g_ratio_phi->GetYaxis()->SetTitle("Data/MC");
    g_ratio_phi->GetYaxis()->SetNdivisions(5);
    g_ratio_phi->GetYaxis()->SetTitleSize(15);
    g_ratio_phi->GetYaxis()->SetTitleFont(43);
    g_ratio_phi->GetYaxis()->SetTitleOffset(1.55);
    g_ratio_phi->GetYaxis()->SetLabelFont(43); 
    g_ratio_phi->GetYaxis()->SetLabelSize(15);
    g_ratio_phi->GetXaxis()->SetTitleSize(15);
    g_ratio_phi->GetXaxis()->SetTitleFont(43);
    g_ratio_phi->GetXaxis()->SetTitleOffset(5.);
    g_ratio_phi->GetXaxis()->SetLabelFont(43);
    g_ratio_phi->GetXaxis()->SetLabelSize(15);
    TLine *oneline_phi = new TLine(g_ratio_phi->GetXaxis()->GetXmin(),1,g_ratio_phi->GetXaxis()->GetXmax(),1);
    oneline_phi->SetLineColor(mit_gray);
    oneline_phi->SetLineWidth(1);
    oneline_phi->SetLineStyle(3);
    g_ratio_phi->Draw("ap");
    oneline_phi->Draw("SAME");
    g_ratio_phi->Draw("p same");
    canvas_phi->Print((data_dir+"eff_ratio_phi.png").c_str());
  }
  if(do_npv) {
    TGraphAsymmErrors *g_data_npv = (TGraphAsymmErrors*) f_data->Get("grEffNPV");
    TGraphAsymmErrors *g_mc_npv = (TGraphAsymmErrors*) f_mc->Get("grEffNPV");
    TGraphAsymmErrors *g_ratio_npv = (TGraphAsymmErrors*) g_data_npv->Clone("g_ratio_npv");
    g_data_npv->SetTitle("");
    g_mc_npv  ->SetTitle("");
    g_ratio_npv->SetTitle(""); 
    unsigned int np_npv = g_data_npv->GetN();
    for(unsigned int ip = 0; ip < np_npv; ip++) {
      Double_t x_data,y_data,x_mc,y_mc;
      g_data_npv->GetPoint(ip,x_data,y_data);
      g_mc_npv->GetPoint(ip,x_mc,y_mc);
      g_ratio_npv->SetPoint(ip, x_data, y_data/y_mc);
      g_ratio_npv->SetPointError(ip,
        g_data_npv->GetErrorXlow(ip),
        g_data_npv->GetErrorXhigh(ip),
        sqrt(pow(g_data_npv->GetErrorYlow(ip)/y_mc,2) + pow(g_mc_npv->GetErrorYhigh(ip)*y_data/(y_mc*y_mc),2)),
        sqrt(pow(g_data_npv->GetErrorYhigh(ip)/y_mc,2) + pow(g_mc_npv->GetErrorYlow(ip)*y_data/(y_mc*y_mc),2))
      );
    }
    logx=false;
    TCanvas *canvas_npv = new TCanvas("canvas_npv", "canvas_npv", 640,480);
    canvas_npv->SetMargin(0,0,0,0);
    TPad *pad1_npv = new TPad("pad1_npv", "pad1_npv", 0, .3, 1, 1);
    pad1_npv->SetGrid(0,1);
    pad1_npv->SetMargin(0.1,0.04,0,.1);
    pad1_npv->Draw();
    pad1_npv->cd();
    if(logx) pad1_npv->SetLogx();
    g_mc_npv->GetYaxis()->SetTitleSize(15);
    g_mc_npv->GetYaxis()->SetTitleFont(43);
    g_mc_npv->GetYaxis()->SetTitleOffset(1.55);
    g_mc_npv->SetMarkerColor(mit_red);
    g_mc_npv->SetLineColor(mit_red);
    g_data_npv->GetYaxis()->SetTitle("#varepsilon");
    g_data_npv->Draw("ap");
    g_data_npv->SetMaximum(1.1);
    g_data_npv->SetMinimum(0.15);
    g_mc_npv->Draw("p same");
    TPad *pad2_npv = new TPad("pad2_npv", "pad2_npv", 0, 0.05, 1, 0.3);
    canvas_npv->cd();
    pad2_npv->SetMargin(0.1,0.04,0.3,0);
    pad2_npv->Draw();
    pad2_npv->cd();
    if(logx) pad2_npv->SetLogx();
    g_ratio_npv->SetMaximum(1.2);
    g_ratio_npv->SetMinimum(0.8);
    if(logx) g_ratio_npv->GetXaxis()->SetMoreLogLabels();
    g_ratio_npv->GetXaxis()->SetTitle("Number of primary vertices");
    g_ratio_npv->GetYaxis()->SetTitle("Data/MC");
    g_ratio_npv->GetYaxis()->SetNdivisions(5);
    g_ratio_npv->GetYaxis()->SetTitleSize(15);
    g_ratio_npv->GetYaxis()->SetTitleFont(43);
    g_ratio_npv->GetYaxis()->SetTitleOffset(1.55);
    g_ratio_npv->GetYaxis()->SetLabelFont(43); 
    g_ratio_npv->GetYaxis()->SetLabelSize(15);
    g_ratio_npv->GetXaxis()->SetTitleSize(15);
    g_ratio_npv->GetXaxis()->SetTitleFont(43);
    g_ratio_npv->GetXaxis()->SetTitleOffset(5.);
    g_ratio_npv->GetXaxis()->SetLabelFont(43);
    g_ratio_npv->GetXaxis()->SetLabelSize(15);
    TLine *oneline_npv = new TLine(g_ratio_npv->GetXaxis()->GetXmin(),1,g_ratio_npv->GetXaxis()->GetXmax(),1);
    oneline_npv->SetLineColor(mit_gray);
    oneline_npv->SetLineWidth(1);
    oneline_npv->SetLineStyle(3);
    g_ratio_npv->Draw("ap");
    oneline_npv->Draw("SAME");
    g_ratio_npv->Draw("p same");
    canvas_npv->Print((data_dir+"eff_ratio_npv.png").c_str());
  }
  TH2D *h_data_etapt_eff, *h_data_etapt_error_hi, *h_data_etapt_error_lo,
       *h_mc_etapt_eff, *h_mc_etapt_error_hi, *h_mc_etapt_error_lo;
  if(do_etapt) {
    if(!do_pt) g_data_pt  = (TGraphAsymmErrors*) f_data->Get("grEffPt");
    h_data_etapt_eff      = (TH2D*) f_data->Get("hEffEtaPt");
    h_data_etapt_error_hi = (TH2D*) f_data->Get("hErrhEtaPt");
    h_data_etapt_error_lo = (TH2D*) f_data->Get("hErrlEtaPt");
    h_mc_etapt_eff      = (TH2D*) f_mc->Get("hEffEtaPt");
    h_mc_etapt_error_hi = (TH2D*) f_mc->Get("hErrhEtaPt");
    h_mc_etapt_error_lo = (TH2D*) f_mc->Get("hErrlEtaPt");
    unsigned int nbins_eta = h_data_etapt_eff->GetNbinsX();
    unsigned int nbins_pt = h_data_etapt_eff->GetNbinsY();
    for(unsigned int nbin_eta = 1; nbin_eta<=nbins_eta; nbin_eta++) {
      TGraphAsymmErrors *g_data_pt_etaslice = new TGraphAsymmErrors(nbins_pt);
      TGraphAsymmErrors *g_mc_pt_etaslice = new TGraphAsymmErrors(nbins_pt);
      TGraphAsymmErrors *g_ratio_pt_etaslice = new TGraphAsymmErrors(nbins_pt);
      for(unsigned int nbin_pt=1; nbin_pt<=nbins_pt; nbin_pt++) {
        unsigned int nbin_2d = h_data_etapt_eff->GetBin(nbin_eta, nbin_pt);
        double data_eff      = h_data_etapt_eff->GetBinContent(nbin_2d),
               data_error_hi = h_data_etapt_error_hi->GetBinContent(nbin_2d),
               data_error_lo = h_data_etapt_error_lo->GetBinContent(nbin_2d),
               mc_eff      = h_mc_etapt_eff->GetBinContent(nbin_2d),
               mc_error_hi = h_mc_etapt_error_hi->GetBinContent(nbin_2d),
               mc_error_lo = h_mc_etapt_error_lo->GetBinContent(nbin_2d);
        double x_data, y_data;
        g_data_pt->GetPoint(nbin_pt-1, x_data, y_data);
        g_data_pt_etaslice->SetPoint(nbin_pt-1, x_data, data_eff);
        g_data_pt_etaslice->SetPointError(nbin_pt-1,
          g_data_pt->GetErrorXlow(nbin_pt-1),
          g_data_pt->GetErrorXhigh(nbin_pt-1),
          data_error_lo,
          data_error_hi
        );
        g_mc_pt_etaslice->SetPoint(nbin_pt-1, x_data, mc_eff);
        g_mc_pt_etaslice->SetPointError(nbin_pt-1,
          g_data_pt->GetErrorXlow(nbin_pt-1),
          g_data_pt->GetErrorXhigh(nbin_pt-1),
          mc_error_lo,
          mc_error_hi
        );
        g_ratio_pt_etaslice->SetPoint(nbin_pt-1, x_data, data_eff/mc_eff);
        g_ratio_pt_etaslice->SetPointError(nbin_pt-1,
          g_data_pt->GetErrorXlow(nbin_pt-1),
          g_data_pt->GetErrorXhigh(nbin_pt-1),
          sqrt(pow(data_error_hi/mc_eff,2) + pow(mc_error_lo*data_eff/(mc_eff*mc_eff),2)),
          sqrt(pow(data_error_lo/mc_eff,2) + pow(mc_error_hi*data_eff/(mc_eff*mc_eff),2))
        );
      }
      char title[512];
      sprintf(title, "#varepsilon as p_{T} for |#eta| [%4.3f, %4.3f]", h_data_etapt_eff->GetXaxis()->GetBinLowEdge(nbin_eta), h_data_etapt_eff->GetXaxis()->GetBinUpEdge(nbin_eta));
      g_data_pt_etaslice->SetTitle(title);
      logx=false;
      TCanvas *canvas_pt_etaslice = new TCanvas("canvas_pt_etaslice", "canvas_pt_etaslice", 640,480);
      canvas_pt_etaslice->SetMargin(0,0,0,0);
      TPad *pad1_pt_etaslice = new TPad("pad1_pt_etaslice", "pad1_pt_etaslice", 0, .3, 1, 1);
      pad1_pt_etaslice->SetGrid(0,1);
      pad1_pt_etaslice->SetMargin(0.1,0.04,0,.1);
      pad1_pt_etaslice->Draw();
      pad1_pt_etaslice->cd();
      if(logx) pad1_pt_etaslice->SetLogx();
      g_mc_pt_etaslice->GetYaxis()->SetTitleSize(15);
      g_mc_pt_etaslice->GetYaxis()->SetTitleFont(43);
      g_mc_pt_etaslice->GetYaxis()->SetTitleOffset(1.55);
      g_mc_pt_etaslice->SetMarkerColor(mit_red);
      g_mc_pt_etaslice->SetMarkerStyle(20);
      g_mc_pt_etaslice->SetLineColor(mit_red);
      g_data_pt_etaslice->GetYaxis()->SetTitle("#varepsilon");
      g_data_pt_etaslice->SetMarkerStyle(20);
      g_data_pt_etaslice->Draw("ap");
      g_data_pt_etaslice->SetMaximum(1.1);
      g_data_pt_etaslice->SetMinimum(0.15);
      g_mc_pt_etaslice->Draw("p same");
      TPad *pad2_pt_etaslice = new TPad("pad2_pt_etaslice", "pad2_pt_etaslice", 0, 0.05, 1, 0.3);
      canvas_pt_etaslice->cd();
      pad2_pt_etaslice->SetMargin(0.1,0.04,0.3,0);
      pad2_pt_etaslice->Draw();
      pad2_pt_etaslice->cd();
      if(logx) pad2_pt_etaslice->SetLogx();
      g_ratio_pt_etaslice->SetMaximum(1.2);
      g_ratio_pt_etaslice->SetMinimum(0.8);
      if(logx) g_ratio_pt_etaslice->GetXaxis()->SetMoreLogLabels();
      g_ratio_pt_etaslice->SetMarkerStyle(20);
      g_ratio_pt_etaslice->SetTitle("");
      g_ratio_pt_etaslice->GetXaxis()->SetTitle("p_{T} [GeV]");
      g_ratio_pt_etaslice->GetYaxis()->SetTitle("Data/MC");
      g_ratio_pt_etaslice->GetYaxis()->SetNdivisions(5);
      g_ratio_pt_etaslice->GetYaxis()->SetTitleSize(15);
      g_ratio_pt_etaslice->GetYaxis()->SetTitleFont(43);
      g_ratio_pt_etaslice->GetYaxis()->SetTitleOffset(1.55);
      g_ratio_pt_etaslice->GetYaxis()->SetLabelFont(43); 
      g_ratio_pt_etaslice->GetYaxis()->SetLabelSize(15);
      g_ratio_pt_etaslice->GetXaxis()->SetTitleSize(15);
      g_ratio_pt_etaslice->GetXaxis()->SetTitleFont(43);
      g_ratio_pt_etaslice->GetXaxis()->SetTitleOffset(5.);
      g_ratio_pt_etaslice->GetXaxis()->SetLabelFont(43);
      g_ratio_pt_etaslice->GetXaxis()->SetLabelSize(15);
      TLine *oneline_pt_etaslice = new TLine(g_ratio_pt_etaslice->GetXaxis()->GetXmin(),1,g_ratio_pt_etaslice->GetXaxis()->GetXmax(),1);
      oneline_pt_etaslice->SetLineColor(mit_gray);
      oneline_pt_etaslice->SetLineWidth(1);
      oneline_pt_etaslice->SetLineStyle(3);
      g_ratio_pt_etaslice->Draw("ap");
      oneline_pt_etaslice->Draw("SAME");
      g_ratio_pt_etaslice->Draw("p same");
      char filename[512];
      sprintf(filename, "%seff_ratio_pt_etaslice_%d.png", data_dir.c_str(), nbin_eta-1);
      canvas_pt_etaslice->Print(filename);
      delete canvas_pt_etaslice;
    }
  }

}

void make_some_ratio_plots_80x() {
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToTight_electronTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToMedium_electronTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToLoose_electronTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToVeto_electronTnP/");
  
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleElectron_BaselineToTight_signEta_electronTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToTight_signEta_electronTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleElectron_BaselineToMedium_signEta_electronTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToMedium_signEta_electronTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleElectron_BaselineToLoose_signEta_electronTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToLoose_signEta_electronTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleElectron_BaselineToVeto_signEta_electronTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToVeto_signEta_electronTnP/",false,true,false,false,false);

  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleMuon_BaselineToTight_coarseEta_muonTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToTight_coarseEta_muonTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleMuon_BaselineToMedium_coarseEta_muonTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToMedium_coarseEta_muonTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleMuon_BaselineToLoose_coarseEta_muonTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToLoose_coarseEta_muonTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleMuon_BaselineToTight_fineEta_muonTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToTight_fineEta_muonTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleMuon_BaselineToMedium_fineEta_muonTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToMedium_fineEta_muonTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/SingleMuon_BaselineToLoose_fineEta_muonTnP/", "~dhsu/TagAndProbe/2016-06-08_80x_500pb/template_erfcexp/DYJetsToLL_BaselineToLoose_fineEta_muonTnP/",false,true,false,false,false);

}
void make_some_ratio_plots_76x() {
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_BaselineToTight_electronTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_BaselineToMedium_electronTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_BaselineToLoose_electronTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_BaselineToVeto_electronTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_signEta_electronTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_BaselineToTight_signEta_electronTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_signEta_electronTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_BaselineToMedium_signEta_electronTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_signEta_electronTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_BaselineToLoose_signEta_electronTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_signEta_electronTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_BaselineToVeto_signEta_electronTnP/",false,true,false,false,false);

  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_muonTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/DYJetsToLL_BaselineToTight_muonTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_muonTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/DYJetsToLL_BaselineToMedium_muonTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_muonTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/DYJetsToLL_BaselineToLoose_muonTnP/");
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_muons_fineEta/template_erfcexp/SingleMuon_BaselineToTight_muonTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_muons_fineEta/template_erfcexp/DYJetsToLL_BaselineToTight_muonTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_muons_fineEta/template_erfcexp/SingleMuon_BaselineToMedium_muonTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_muons_fineEta/template_erfcexp/DYJetsToLL_BaselineToMedium_muonTnP/",false,true,false,false,false);
  plot_sf_1d("~dhsu/TagAndProbe/2016-05-02_76x_muons_fineEta/template_erfcexp/SingleMuon_BaselineToLoose_muonTnP/", "~dhsu/TagAndProbe/2016-05-02_76x_muons_fineEta/template_erfcexp/DYJetsToLL_BaselineToLoose_muonTnP/",false,true,false,false,false);

}
