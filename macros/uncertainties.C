#include <TROOT.h>
#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <vector>
#include <fstream>
#include <sstream>
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
using namespace std;

Float_t ele_pt_bins[] = {10,20,30,40,50,200};
Float_t ele_eta_bins[] = {0., 0.8, 1.4442, 1.566, 2., 2.5};
Int_t n_ele_pt_bins=5;
Int_t n_ele_eta_bins=5;
Float_t mu_pt_bins[] = {10,20,30,40,50,70,100,8000};
Float_t mu_eta_bins[] = {0, 1.479, 2.4};
Int_t n_mu_pt_bins=7;
Int_t n_mu_eta_bins=2;

Int_t mit_red  = 1861; 
Int_t mit_gray = 1862; 

void data_style (TH1D *histo) {
  histo->SetLineColor(1);
  histo->SetMarkerColor(1);
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(0.8);
}
void mc_style (TH1D *histo) {
  //histo->SetMinimum(.1);
  //histo->SetMaximum(pow(10,5));
  histo->SetFillColor(mit_red);
  histo->SetFillStyle(1001);
  histo->SetLineColor(1);
}

TCanvas *ratio_plot(
  string basename,
  string canvas_title,
  string xlabel,
  TH1D *histo_data,
  TH1D *histo_mc,
  bool logx = false,
  bool logy = true,
  bool restyle_histos = true
) {
  TCanvas *canvas = new TCanvas(("c_"+basename).c_str(), canvas_title.c_str(), 800,600);
  canvas->SetMargin(0,0,0,0);
  TPad *pad1 = new TPad(("pad1_"+basename).c_str(), "pad1", 0, .3, 1, 1);
  pad1->SetGrid(0,1);
  pad1->SetMargin(0.1,0.04,0.07,.1);
  pad1->Draw();
  pad1->cd();
  if(logx) pad1->SetLogx();
  if(logy) pad1->SetLogy();
  if(restyle_histos) mc_style(histo_mc);
  if(restyle_histos) data_style(histo_data);
  histo_mc->SetTitle(canvas_title.c_str());
  if(restyle_histos) histo_mc->Draw("B HIST");
  else histo_mc->Draw("P0 E1");
  if(logx) histo_mc->GetXaxis()->SetMoreLogLabels();
  histo_mc->GetYaxis()->SetTitleSize(15);
  histo_mc->GetYaxis()->SetTitleFont(43);
  histo_mc->GetYaxis()->SetTitleOffset(1.55);
  //histo_mc->GetYaxis()->SetTitle("Events");
  if(restyle_histos) histo_data->Draw("P E0 X0 SAME");
  else histo_data->Draw("P0 E1 SAME");

  TLegend *legend = new TLegend(0.75, 0.7, 0.92, 0.85);
  legend->AddEntry(histo_data, "Data", "lp");
  if(restyle_histos) legend->AddEntry(histo_mc, "DY MC", "f");
  else legend->AddEntry(histo_mc, "DY MC", "lp");
  legend->SetFillColor(0);
  legend->Draw("SAME");
  TPad *pad2 = new TPad(("pad2_"+basename).c_str(), "pad2", 0, 0.05, 1, 0.3);
  canvas->cd();
  pad2->SetMargin(0.1,0.04,0.3,0.04);
  pad2->Draw();
  pad2->cd();
  if(logx) pad2->SetLogx();
  TH1D *histo_ratio = (TH1D*)histo_data->Clone();
  histo_ratio->SetTitle("");;
  histo_ratio->Divide( histo_mc);
  data_style(histo_ratio);
  histo_ratio->SetMaximum(1.4);
  histo_ratio->SetMinimum(0.6);
  histo_ratio->Draw("P E0 X0");
  if(logx) histo_ratio->GetXaxis()->SetMoreLogLabels();
  histo_ratio->GetXaxis()->SetTitle(xlabel.c_str());
  histo_ratio->GetYaxis()->SetTitle("Data/MC");
  histo_ratio->GetYaxis()->SetNdivisions(5);
  histo_ratio->GetYaxis()->SetTitleSize(15);
  histo_ratio->GetYaxis()->SetTitleFont(43);
  histo_ratio->GetYaxis()->SetTitleOffset(1.55);
  histo_ratio->GetYaxis()->SetLabelFont(43); 
  histo_ratio->GetYaxis()->SetLabelSize(15);
  histo_ratio->GetXaxis()->SetTitleSize(15);
  histo_ratio->GetXaxis()->SetTitleFont(43);
  histo_ratio->GetXaxis()->SetTitleOffset(4.);
  histo_ratio->GetXaxis()->SetLabelFont(43);
  histo_ratio->GetXaxis()->SetLabelSize(15);
  double xlo  = histo_data->GetXaxis()->GetBinLowEdge(1);
  double xhi  = histo_data->GetXaxis()->GetBinUpEdge(histo_data->GetNbinsX());
  
  TLine *oneline = new TLine(xlo,1,xhi,1);
  oneline->SetLineColor(1);
  oneline->SetLineWidth(1);
  oneline->Draw("SAME");
  printf("%s : data integral %f, mc integral %f \n", basename.c_str(), histo_data->Integral(), histo_mc->Integral());
  return canvas;

}





void uncertainties(
  string plots_dir,
  string root_dir,
  bool draw = false
) {
  gStyle->SetOptStat(0); 
  //open files
  TFile *f_mu_sf_BWCBPlusVoigt_erfcexp       = TFile::Open("~/leptonScaleFactors/root/16-02-2016/BWCBPlusVoigt_erfcexp/scalefactors_mu.root","READ");
  TFile *f_mu_sf_template_erfcexp = TFile::Open("~/leptonScaleFactors/root/16-02-2016/template_erfcexp/scalefactors_mu.root","READ");
  TFile *f_mu_sf_template_exp     = TFile::Open("~/leptonScaleFactors/root/16-02-2016/template_exp/scalefactors_mu.root","READ");
  TFile *f_ele_sf_BWCBPlusVoigt_erfcexp       = TFile::Open("~/leptonScaleFactors/root/16-02-2016/BWCBPlusVoigt_erfcexp/scalefactors_ele.root","READ");
  TFile *f_ele_sf_template_erfcexp = TFile::Open("~/leptonScaleFactors/root/16-02-2016/template_erfcexp/scalefactors_ele.root","READ");
  TFile *f_ele_sf_template_exp     = TFile::Open("~/leptonScaleFactors/root/16-02-2016/template_exp/scalefactors_ele.root","READ");
  
  TFile *f_ele_eff_data = TFile::Open("~/leptonScaleFactors/root/16-02-2016/template_erfcexp/SingleElectron_efficiencies_electronTnP.root","READ");
  TFile *f_ele_eff_mc = TFile::Open("~/leptonScaleFactors/root/16-02-2016/template_erfcexp/DYJetsToLL_efficiencies_electronTnP.root","READ");
  // systematic uncertainties on factorized electron efficiency from data due to factorization
  TH2D *syst_ele_eff_fact_veto  = (TH2D*) f_ele_eff_data->Get("absdiff_Veto_ele");
  TH2D *syst_ele_eff_fact_tight = (TH2D*) f_ele_eff_data->Get("absdiff_Tight_ele");
  // electron efficiencies in MC
  TH2D *ele_eff_MC_veto  = (TH2D*) f_ele_eff_mc->Get("factEff_Veto_ele");
  TH2D *ele_eff_MC_tight = (TH2D*) f_ele_eff_mc->Get("factEff_Tight_ele");
  
  // divide the two to get the systematic uncertainty on the electron scale factors due to factorization
  TH2D *syst_ele_sf_fact_veto  = new TH2D("syst_ele_sf_fact_veto",  "syst_ele_sf_fact_veto",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *syst_ele_sf_fact_tight = new TH2D("syst_ele_sf_fact_tight", "syst_ele_sf_fact_Right", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  *syst_ele_sf_fact_veto = (*syst_ele_eff_fact_veto) / (*ele_eff_MC_veto);
  *syst_ele_sf_fact_tight = (*syst_ele_eff_fact_tight) / (*ele_eff_MC_tight);

  // systematic uncertainty from signal model choice
  TH2D *syst_ele_sf_signal_veto  = (TH2D*) f_ele_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_factorized_scalefactors_Veto_ele");
  TH2D *syst_ele_sf_signal_tight = (TH2D*) f_ele_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_factorized_scalefactors_Tight_ele");
  TH2D *syst_mu_sf_signal_loose  = (TH2D*) f_mu_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Loose_mu");
  TH2D *syst_mu_sf_signal_tight = (TH2D*) f_mu_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Tight_mu");
  
  // systematic uncertainty from background model choice
  TH2D *syst_ele_sf_background_veto  = (TH2D*) f_ele_sf_template_exp->Get("absDiff_factorized_scalefactors_Veto_ele");
  TH2D *syst_ele_sf_background_tight = (TH2D*) f_ele_sf_template_exp->Get("absDiff_factorized_scalefactors_Tight_ele");
  TH2D *syst_mu_sf_background_loose = (TH2D*) f_mu_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Loose_mu");
  TH2D *syst_mu_sf_background_tight = (TH2D*) f_mu_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Tight_mu");

  // combined systematic uncertainty
  TH2D *syst_ele_sf_combined_veto   = new TH2D("syst_ele_sf_combined_veto", "Combined syst. unc. on electron veto scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *syst_ele_sf_combined_tight  = new TH2D("syst_ele_sf_combined_tight", "Combined syst. unc. on electron tight scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *syst_mu_sf_combined_loose  = new TH2D("syst_mu_sf_combined_loose", "Combined syst. unc. on muon loose scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  TH2D *syst_mu_sf_combined_tight  = new TH2D("syst_mu_sf_combined_tight", "Combined syst. unc. on muon tight scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  
  // statistical uncertainty
  TH2D *ele_sf_veto        = (TH2D*) f_ele_sf_template_erfcexp ->Get("factorized_scalefactors_Veto_ele"); 
  TH2D *ele_sf_tight       = (TH2D*) f_ele_sf_template_erfcexp ->Get("factorized_scalefactors_Tight_ele"); 
  TH2D *mu_sf_loose        = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Loose_mu");
  TH2D *mu_sf_tight        = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Tight_mu");
  TH2D *stat_ele_sf_veto   = new TH2D("stat_ele_sf_veto",  "Stat. unc. on electron veto scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *stat_ele_sf_tight  = new TH2D("stat_ele_sf_tight", "Stat. unc. on electron tight scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *stat_mu_sf_loose   = new TH2D("stat_mu_sf_loose",  "Stat. unc. on muon loose scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  TH2D *stat_mu_sf_tight   = new TH2D("stat_mu_sf_tight",  "Stat. unc. on muon tight scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  
  for(int i_eta = 1; i_eta <= n_ele_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_ele_pt_bins; i_pt++) {
    int n = syst_ele_sf_combined_veto->GetBin(i_eta, i_pt);
    // assign combined electron syst. unc.
    syst_ele_sf_combined_veto->SetBinContent(n, sqrt(
      pow(syst_ele_sf_fact_veto       ->GetBinContent(n), 2) + 
      pow(syst_ele_sf_signal_veto     ->GetBinContent(n), 2) + 
      pow(syst_ele_sf_background_veto ->GetBinContent(n), 2) 
    )); 
    syst_ele_sf_combined_tight->SetBinContent(n, sqrt(
      pow(syst_ele_sf_fact_tight       ->GetBinContent(n), 2) + 
      pow(syst_ele_sf_signal_tight     ->GetBinContent(n), 2) + 
      pow(syst_ele_sf_background_tight ->GetBinContent(n), 2) 
    )); 
    // pull out electron stat. unc
    stat_ele_sf_veto->SetBinContent(n, ele_sf_veto->GetBinError(n));
    stat_ele_sf_tight->SetBinContent(n, ele_sf_tight->GetBinError(n));
  }}
  for(int i_eta = 1; i_eta <= n_mu_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_mu_pt_bins; i_pt++) {
    int n = syst_mu_sf_signal_loose->GetBin(i_eta, i_pt);
    // assign combined muon syst. unc.
    syst_mu_sf_combined_loose->SetBinContent(n, sqrt(
      pow(syst_mu_sf_signal_loose     ->GetBinContent(n), 2) +
      pow(syst_mu_sf_background_loose ->GetBinContent(n), 2)     
    )); 
    syst_mu_sf_combined_tight->SetBinContent(n, sqrt(
      pow(syst_mu_sf_signal_tight     ->GetBinContent(n), 2) + 
      pow(syst_mu_sf_background_tight ->GetBinContent(n), 2) 
    )); 
    // pull out electron stat. unc
    stat_mu_sf_loose->SetBinContent(n, mu_sf_loose->GetBinError(n));
    stat_mu_sf_tight->SetBinContent(n, mu_sf_tight->GetBinError(n));
  }}
  
  if(draw) {
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.3f");
    mitPalette();
    TPaletteAxis *palette_axis;

    // ele veto plots
    TCanvas *c_syst_ele_sf_fact_veto = new TCanvas("c_syst_ele_sf_fact_veto", "Syst. unc. on ele. veto SF (fact.)", 400, 600);
    syst_ele_sf_fact_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_fact_veto->SetTitle("Syst. unc. on ele. veto SF from factorization");
    syst_ele_sf_fact_veto->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_fact_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_fact_veto->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_fact_veto->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_fact_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_fact_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_fact_veto->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_fact_veto->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_fact_veto->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_fact_veto->SetMinimum(0);
    syst_ele_sf_fact_veto->SetMaximum(0.2);
    syst_ele_sf_fact_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_fact_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_fact_veto->Update();
    c_syst_ele_sf_fact_veto->Print((plots_dir+"syst_ele_sf_fact_veto.png").c_str());
    
    TCanvas *c_syst_ele_sf_signal_veto = new TCanvas("c_syst_ele_sf_signal_veto", "Syst. unc. on ele. veto SF (signal)", 400, 600);
    syst_ele_sf_signal_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_signal_veto->SetTitle("Syst. unc. on ele. veto SF from signal choice");
    syst_ele_sf_signal_veto->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_signal_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_veto->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_veto->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_signal_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_veto->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_veto->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_veto->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_signal_veto->SetMinimum(0);
    syst_ele_sf_signal_veto->SetMaximum(0.2);
    syst_ele_sf_signal_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_signal_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_signal_veto->Update();
    c_syst_ele_sf_signal_veto->Print((plots_dir+"syst_ele_sf_signal_veto.png").c_str());
    
    TCanvas *c_syst_ele_sf_background_veto = new TCanvas("c_syst_ele_sf_background_veto", "Syst. unc. on ele. veto SF (b.g.)", 400, 600);
    syst_ele_sf_background_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_background_veto->SetTitle("Syst. unc. on ele. veto SF from background choice");
    syst_ele_sf_background_veto->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_background_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_veto->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_veto->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_background_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_veto->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_veto->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_veto->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_background_veto->SetMinimum(0);
    syst_ele_sf_background_veto->SetMaximum(0.2);
    syst_ele_sf_background_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_background_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_background_veto->Update();
    c_syst_ele_sf_background_veto->Print((plots_dir+"syst_ele_sf_background_veto.png").c_str());
    
    TCanvas *c_syst_ele_sf_combined_veto = new TCanvas("c_syst_ele_sf_combined_veto", "Syst. unc. on ele. veto SF (combined)", 400, 600);
    syst_ele_sf_combined_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_combined_veto->SetTitle("Syst. unc. on ele. veto SF (combined)");
    syst_ele_sf_combined_veto->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_combined_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_veto->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_veto->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_combined_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_veto->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_veto->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_veto->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_combined_veto->SetMinimum(0);
    syst_ele_sf_combined_veto->SetMaximum(0.2);
    syst_ele_sf_combined_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_combined_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_combined_veto->Update();
    c_syst_ele_sf_combined_veto->Print((plots_dir+"syst_ele_sf_combined_veto.png").c_str());
    
    TCanvas *c_stat_ele_sf_veto = new TCanvas("c_stat_ele_sf_veto", "Stat. unc. on ele. veto SF", 400, 600);
    stat_ele_sf_veto->Draw("TEXT COLZ");
    gPad->Update();
    stat_ele_sf_veto->SetTitle("Stat. unc. on ele. veto SF");
    stat_ele_sf_veto->GetXaxis()->SetTitle("| #eta |");
    stat_ele_sf_veto->GetXaxis()->SetTitleOffset(0.9);
    stat_ele_sf_veto->GetXaxis()->SetTitleSize(0.04);
    stat_ele_sf_veto->GetXaxis()->SetLabelSize(0.02);
    stat_ele_sf_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_ele_sf_veto->GetYaxis()->SetTitleOffset(0.9);
    stat_ele_sf_veto->GetYaxis()->SetTitleSize(0.04);
    stat_ele_sf_veto->GetYaxis()->SetLabelSize(0.02);
    stat_ele_sf_veto->GetYaxis()->SetRangeUser(10,100);
    stat_ele_sf_veto->SetMinimum(0);
    stat_ele_sf_veto->SetMaximum(0.2);
    stat_ele_sf_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_ele_sf_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_ele_sf_veto->Update();
    c_stat_ele_sf_veto->Print((plots_dir+"stat_ele_sf_veto.png").c_str());

    // ele tight plots
    TCanvas *c_syst_ele_sf_fact_tight = new TCanvas("c_syst_ele_sf_fact_tight", "Syst. unc. on ele. tight SF (fact.)", 400, 600);
    syst_ele_sf_fact_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_fact_tight->SetTitle("Syst. unc. on ele. tight SF from factorization");
    syst_ele_sf_fact_tight->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_fact_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_fact_tight->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_fact_tight->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_fact_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_fact_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_fact_tight->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_fact_tight->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_fact_tight->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_fact_tight->SetMinimum(0);
    syst_ele_sf_fact_tight->SetMaximum(0.2);
    syst_ele_sf_fact_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_fact_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_fact_tight->Update();
    c_syst_ele_sf_fact_tight->Print((plots_dir+"syst_ele_sf_fact_tight.png").c_str());
    
    TCanvas *c_syst_ele_sf_signal_tight = new TCanvas("c_syst_ele_sf_signal_tight", "Syst. unc. on ele. tight SF (signal)", 400, 600);
    syst_ele_sf_signal_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_signal_tight->SetTitle("Syst. unc. on ele. tight SF from signal choice");
    syst_ele_sf_signal_tight->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_signal_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_tight->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_tight->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_signal_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_tight->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_tight->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_tight->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_signal_tight->SetMinimum(0);
    syst_ele_sf_signal_tight->SetMaximum(0.2);
    syst_ele_sf_signal_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_signal_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_signal_tight->Update();
    c_syst_ele_sf_signal_tight->Print((plots_dir+"syst_ele_sf_signal_tight.png").c_str());
    
    TCanvas *c_syst_ele_sf_background_tight = new TCanvas("c_syst_ele_sf_background_tight", "Syst. unc. on ele. tight SF (b.g.)", 400, 600);
    syst_ele_sf_background_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_background_tight->SetTitle("Syst. unc. on ele. tight SF from background choice");
    syst_ele_sf_background_tight->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_background_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_tight->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_tight->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_background_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_tight->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_tight->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_tight->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_background_tight->SetMinimum(0);
    syst_ele_sf_background_tight->SetMaximum(0.2);
    syst_ele_sf_background_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_background_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_background_tight->Update();
    c_syst_ele_sf_background_tight->Print((plots_dir+"syst_ele_sf_background_tight.png").c_str());
    
    TCanvas *c_syst_ele_sf_combined_tight = new TCanvas("c_syst_ele_sf_combined_tight", "Syst. unc. on ele. tight SF (combined)", 400, 600);
    syst_ele_sf_combined_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_combined_tight->SetTitle("Syst. unc. on ele. tight SF (combined)");
    syst_ele_sf_combined_tight->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_combined_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_tight->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_tight->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_combined_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_tight->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_tight->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_tight->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_combined_tight->SetMinimum(0);
    syst_ele_sf_combined_tight->SetMaximum(0.2);
    syst_ele_sf_combined_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_combined_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_combined_tight->Update();
    c_syst_ele_sf_combined_tight->Print((plots_dir+"syst_ele_sf_combined_tight.png").c_str());
    
    TCanvas *c_stat_ele_sf_tight = new TCanvas("c_stat_ele_sf_tight", "Stat. unc. on ele. tight SF", 400, 600);
    stat_ele_sf_tight->Draw("TEXT COLZ");
    gPad->Update();
    stat_ele_sf_tight->SetTitle("Stat. unc. on ele. tight SF");
    stat_ele_sf_tight->GetXaxis()->SetTitle("| #eta |");
    stat_ele_sf_tight->GetXaxis()->SetTitleOffset(0.9);
    stat_ele_sf_tight->GetXaxis()->SetTitleSize(0.04);
    stat_ele_sf_tight->GetXaxis()->SetLabelSize(0.02);
    stat_ele_sf_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_ele_sf_tight->GetYaxis()->SetTitleOffset(0.9);
    stat_ele_sf_tight->GetYaxis()->SetTitleSize(0.04);
    stat_ele_sf_tight->GetYaxis()->SetLabelSize(0.02);
    stat_ele_sf_tight->GetYaxis()->SetRangeUser(10,100);
    stat_ele_sf_tight->SetMinimum(0);
    stat_ele_sf_tight->SetMaximum(0.2);
    stat_ele_sf_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_ele_sf_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_ele_sf_tight->Update();
    c_stat_ele_sf_tight->Print((plots_dir+"stat_ele_sf_tight.png").c_str());

    // mu loose plots
    TCanvas *c_syst_mu_sf_signal_loose = new TCanvas("c_syst_mu_sf_signal_loose", "Syst. unc. on mu. loose SF (signal)", 400, 600);
    syst_mu_sf_signal_loose->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_signal_loose->SetTitle("Syst. unc. on mu. loose SF from signal choice");
    syst_mu_sf_signal_loose->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_signal_loose->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_loose->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_loose->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_signal_loose->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_loose->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_loose->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_loose->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_signal_loose->SetMinimum(0);
    syst_mu_sf_signal_loose->SetMaximum(0.2);
    syst_mu_sf_signal_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_signal_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_signal_loose->Update();
    c_syst_mu_sf_signal_loose->Print((plots_dir+"syst_mu_sf_signal_loose.png").c_str());
    
    TCanvas *c_syst_mu_sf_background_loose = new TCanvas("c_syst_mu_sf_background_loose", "Syst. unc. on mu. loose SF (b.g.)", 400, 600);
    syst_mu_sf_background_loose->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_background_loose->SetTitle("Syst. unc. on mu. loose SF from background choice");
    syst_mu_sf_background_loose->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_background_loose->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_loose->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_loose->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_background_loose->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_loose->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_loose->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_loose->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_background_loose->SetMinimum(0);
    syst_mu_sf_background_loose->SetMaximum(0.2);
    syst_mu_sf_background_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_background_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_background_loose->Update();
    c_syst_mu_sf_background_loose->Print((plots_dir+"syst_mu_sf_background_loose.png").c_str());
    
    TCanvas *c_syst_mu_sf_combined_loose = new TCanvas("c_syst_mu_sf_combined_loose", "Syst. unc. on mu. loose SF (combined)", 400, 600);
    syst_mu_sf_combined_loose->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_combined_loose->SetTitle("Syst. unc. on mu. loose SF (combined)");
    syst_mu_sf_combined_loose->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_combined_loose->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_loose->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_loose->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_combined_loose->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_loose->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_loose->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_loose->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_combined_loose->SetMinimum(0);
    syst_mu_sf_combined_loose->SetMaximum(0.2);
    syst_mu_sf_combined_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_combined_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_combined_loose->Update();
    c_syst_mu_sf_combined_loose->Print((plots_dir+"syst_mu_sf_combined_loose.png").c_str());
    
    TCanvas *c_stat_mu_sf_loose = new TCanvas("c_stat_mu_sf_loose", "Stat. unc. on mu. loose SF", 400, 600);
    stat_mu_sf_loose->Draw("TEXT COLZ");
    gPad->Update();
    stat_mu_sf_loose->SetTitle("Stat. unc. on mu. loose SF");
    stat_mu_sf_loose->GetXaxis()->SetTitle("| #eta |");
    stat_mu_sf_loose->GetXaxis()->SetTitleOffset(0.9);
    stat_mu_sf_loose->GetXaxis()->SetTitleSize(0.04);
    stat_mu_sf_loose->GetXaxis()->SetLabelSize(0.02);
    stat_mu_sf_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_mu_sf_loose->GetYaxis()->SetTitleOffset(0.9);
    stat_mu_sf_loose->GetYaxis()->SetTitleSize(0.04);
    stat_mu_sf_loose->GetYaxis()->SetLabelSize(0.02);
    stat_mu_sf_loose->GetYaxis()->SetRangeUser(10,100);
    stat_mu_sf_loose->SetMinimum(0);
    stat_mu_sf_loose->SetMaximum(0.2);
    stat_mu_sf_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_mu_sf_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_mu_sf_loose->Update();
    c_stat_mu_sf_loose->Print((plots_dir+"stat_mu_sf_loose.png").c_str());

    // mu tight plots
    TCanvas *c_syst_mu_sf_signal_tight = new TCanvas("c_syst_mu_sf_signal_tight", "Syst. unc. on mu. tight SF (signal)", 400, 600);
    syst_mu_sf_signal_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_signal_tight->SetTitle("Syst. unc. on mu. tight SF from signal choice");
    syst_mu_sf_signal_tight->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_signal_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_tight->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_tight->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_signal_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_tight->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_tight->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_tight->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_signal_tight->SetMinimum(0);
    syst_mu_sf_signal_tight->SetMaximum(0.2);
    syst_mu_sf_signal_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_signal_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_signal_tight->Update();
    c_syst_mu_sf_signal_tight->Print((plots_dir+"syst_mu_sf_signal_tight.png").c_str());
    
    TCanvas *c_syst_mu_sf_background_tight = new TCanvas("c_syst_mu_sf_background_tight", "Syst. unc. on mu. tight SF (b.g.)", 400, 600);
    syst_mu_sf_background_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_background_tight->SetTitle("Syst. unc. on mu. tight SF from background choice");
    syst_mu_sf_background_tight->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_background_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_tight->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_tight->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_background_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_tight->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_tight->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_tight->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_background_tight->SetMinimum(0);
    syst_mu_sf_background_tight->SetMaximum(0.2);
    syst_mu_sf_background_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_background_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_background_tight->Update();
    c_syst_mu_sf_background_tight->Print((plots_dir+"syst_mu_sf_background_tight.png").c_str());
    
    TCanvas *c_syst_mu_sf_combined_tight = new TCanvas("c_syst_mu_sf_combined_tight", "Syst. unc. on mu. tight SF (combined)", 400, 600);
    syst_mu_sf_combined_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_combined_tight->SetTitle("Syst. unc. on mu. tight SF (combined)");
    syst_mu_sf_combined_tight->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_combined_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_tight->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_tight->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_combined_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_tight->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_tight->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_tight->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_combined_tight->SetMinimum(0);
    syst_mu_sf_combined_tight->SetMaximum(0.2);
    syst_mu_sf_combined_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_combined_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_combined_tight->Update();
    c_syst_mu_sf_combined_tight->Print((plots_dir+"syst_mu_sf_combined_tight.png").c_str());
    
    TCanvas *c_stat_mu_sf_tight = new TCanvas("c_stat_mu_sf_tight", "Stat. unc. on mu. tight SF", 400, 600);
    stat_mu_sf_tight->Draw("TEXT COLZ");
    gPad->Update();
    stat_mu_sf_tight->SetTitle("Stat. unc. on mu. tight SF");
    stat_mu_sf_tight->GetXaxis()->SetTitle("| #eta |");
    stat_mu_sf_tight->GetXaxis()->SetTitleOffset(0.9);
    stat_mu_sf_tight->GetXaxis()->SetTitleSize(0.04);
    stat_mu_sf_tight->GetXaxis()->SetLabelSize(0.02);
    stat_mu_sf_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_mu_sf_tight->GetYaxis()->SetTitleOffset(0.9);
    stat_mu_sf_tight->GetYaxis()->SetTitleSize(0.04);
    stat_mu_sf_tight->GetYaxis()->SetLabelSize(0.02);
    stat_mu_sf_tight->GetYaxis()->SetRangeUser(10,100);
    stat_mu_sf_tight->SetMinimum(0);
    stat_mu_sf_tight->SetMaximum(0.2);
    stat_mu_sf_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_mu_sf_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_mu_sf_tight->Update();
    c_stat_mu_sf_tight->Print((plots_dir+"stat_mu_sf_tight.png").c_str());


  }
  
  TFile *uncertainties_ele_veto = new TFile((root_dir+"combined_uncertainties_ele_veto.root").c_str(),"RECREATE");
  syst_ele_sf_fact_veto       ->Write("syst_ele_sf_fact_veto");
  syst_ele_sf_signal_veto     ->Write("syst_ele_sf_signal_veto");
  syst_ele_sf_background_veto ->Write("syst_ele_sf_background_veto");
  syst_ele_sf_combined_veto   ->Write("syst_ele_sf_combined_veto");
  stat_ele_sf_veto            ->Write("stat_ele_sf_veto");
  TFile *uncertainties_ele_tight = new TFile((root_dir+"combined_uncertainties_ele_tight.root").c_str(),"RECREATE");
  syst_ele_sf_fact_tight       ->Write("syst_ele_sf_fact_tight");
  syst_ele_sf_signal_tight     ->Write("syst_ele_sf_signal_tight");
  syst_ele_sf_background_tight ->Write("syst_ele_sf_background_tight");
  syst_ele_sf_combined_tight   ->Write("syst_ele_sf_combined_tight");
  stat_ele_sf_tight            ->Write("stat_ele_sf_tight");
  TFile *uncertainties_mu_loose = new TFile((root_dir+"combined_uncertainties_mu_loose.root").c_str(),"RECREATE");
  syst_mu_sf_signal_loose     ->Write("syst_mu_sf_signal_loose");
  syst_mu_sf_background_loose ->Write("syst_mu_sf_background_loose");
  syst_mu_sf_combined_loose   ->Write("syst_mu_sf_combined_loose");
  stat_mu_sf_loose            ->Write("stat_mu_sf_loose");
  TFile *uncertainties_mu_tight = new TFile((root_dir+"combined_uncertainties_mu_tight.root").c_str(),"RECREATE");
  syst_mu_sf_signal_tight     ->Write("syst_mu_sf_signal_tight");
  syst_mu_sf_background_tight ->Write("syst_mu_sf_background_tight");
  syst_mu_sf_combined_tight   ->Write("syst_mu_sf_combined_tight");
  stat_mu_sf_tight            ->Write("stat_mu_sf_tight");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void propagate_to_Zpt(
  string plots_dir,
  string root_dir,
  bool verbose=false
) {
  gStyle->SetOptStat(0); 
  int nbins=8;
  int max_pT = 400;
  
  TFile *uncertainties_ele = TFile::Open((root_dir+"scalefactors_ele_76x.root").c_str(),"READ");
  TFile *uncertainties_mu  = TFile::Open((root_dir+"scalefactors_mu_76x.root").c_str(),"READ");
  
  TFile *f_tnp_ele_veto   = TFile::Open("~/leptonScaleFactors/root/DYJetsToLL_NLO_76x_BaselineToVeto_electronTnP.root","READ");
  TFile *f_tnp_ele_tight  = TFile::Open("~/leptonScaleFactors/root/DYJetsToLL_NLO_76x_BaselineToTight_electronTnP.root","READ");
  TFile *f_tnp_mu_loose   = TFile::Open("~/leptonScaleFactors/root/DYJetsToLL_NLO_76x_BaselineToLoose_electronTnP.root","READ");
  TFile *f_tnp_mu_tight   = TFile::Open("~/leptonScaleFactors/root/DYJetsToLL_NLO_76x_BaselineToTight_electronTnP.root","READ");

  if(
    !uncertainties_ele  ||
    !uncertainties_mu   ||
    !f_tnp_ele_veto           ||
    !f_tnp_ele_tight          ||
    !f_tnp_mu_loose           ||
    !f_tnp_mu_tight           
  ) { printf("failed to open file\n"); exit(-1); }
  TH2D *syst_ele_sf_combined_veto         = (TH2D*) uncertainties_ele ->Get("scalefactors_Veto_Electron_syst_error_combined"        );
  TH2D *syst_ele_sf_tag_cuts_veto         = (TH2D*) uncertainties_ele ->Get("scalefactors_Veto_Electron_syst_error_tag_cuts"        );
  TH2D *syst_ele_sf_generator_veto       = (TH2D*) uncertainties_ele ->Get("scalefactors_Veto_Electron_syst_error_generator"       );
  TH2D *syst_ele_sf_background_shape_veto = (TH2D*) uncertainties_ele ->Get("scalefactors_Veto_Electron_syst_error_background_shape");
  TH2D *syst_ele_sf_signal_shape_veto     = (TH2D*) uncertainties_ele ->Get("scalefactors_Veto_Electron_syst_error_signal_shape"    );
  
  TH2D *syst_ele_sf_combined_tight         = (TH2D*) uncertainties_ele ->Get("scalefactors_Tight_Electron_syst_error_combined"        );
  TH2D *syst_ele_sf_tag_cuts_tight         = (TH2D*) uncertainties_ele ->Get("scalefactors_Tight_Electron_syst_error_tag_cuts"        );
  TH2D *syst_ele_sf_generator_tight       = (TH2D*) uncertainties_ele ->Get("scalefactors_Tight_Electron_syst_error_generator"       );
  TH2D *syst_ele_sf_background_shape_tight = (TH2D*) uncertainties_ele ->Get("scalefactors_Tight_Electron_syst_error_background_shape");
  TH2D *syst_ele_sf_signal_shape_tight     = (TH2D*) uncertainties_ele ->Get("scalefactors_Tight_Electron_syst_error_signal_shape"    );
  
  TH2D *syst_mu_sf_combined_loose         = (TH2D*) uncertainties_mu ->Get("scalefactors_Loose_Muon_syst_error_combined"        );
  TH2D *syst_mu_sf_tag_cuts_loose         = (TH2D*) uncertainties_mu ->Get("scalefactors_Loose_Muon_syst_error_tag_cuts"        );
  TH2D *syst_mu_sf_generator_loose       = (TH2D*) uncertainties_mu ->Get("scalefactors_Loose_Muon_syst_error_generator"       );
  TH2D *syst_mu_sf_background_shape_loose = (TH2D*) uncertainties_mu ->Get("scalefactors_Loose_Muon_syst_error_background_shape");
  TH2D *syst_mu_sf_signal_shape_loose     = (TH2D*) uncertainties_mu ->Get("scalefactors_Loose_Muon_syst_error_signal_shape"    );
  
  TH2D *syst_mu_sf_combined_tight         = (TH2D*) uncertainties_mu ->Get("scalefactors_Tight_Muon_syst_error_combined"        );
  TH2D *syst_mu_sf_tag_cuts_tight         = (TH2D*) uncertainties_mu ->Get("scalefactors_Tight_Muon_syst_error_tag_cuts"        );
  TH2D *syst_mu_sf_generator_tight       = (TH2D*) uncertainties_mu ->Get("scalefactors_Tight_Muon_syst_error_generator"       );
  TH2D *syst_mu_sf_background_shape_tight = (TH2D*) uncertainties_mu ->Get("scalefactors_Tight_Muon_syst_error_background_shape");
  TH2D *syst_mu_sf_signal_shape_tight     = (TH2D*) uncertainties_mu ->Get("scalefactors_Tight_Muon_syst_error_signal_shape"    );
  
  // read from tnp skim
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
  
  TH1D *h_Z_pT_ele_veto  = new TH1D("h_Z_pT_ele_veto", "Z p_{T} spectrum for Veto electron pairs",  nbins, 0, max_pT);
  TH1D *h_Z_pT_ele_tight = new TH1D("h_Z_pT_ele_tight", "Z p_{T} spectrum for Tight electron pairs", nbins, 0, max_pT);
  TH1D *h_Z_pT_mu_loose  = new TH1D("h_Z_pT_mu_loose", "Z p_{T} spectrum for Loose muon pairs",  nbins, 0, max_pT);
  TH1D *h_Z_pT_mu_tight  = new TH1D("h_Z_pT_mu_tight", "Z p_{T} spectrum for Tight muon pairs",  nbins, 0, max_pT);
  
  TH1D *h_syst_Z_ele_combined_veto         = new TH1D("h_syst_Z_ele_combined_veto"        ,  "SF syst. unc. as Z p_{T} for Veto electron pairs (combined)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_ele_tag_cuts_veto         = new TH1D("h_syst_Z_ele_tag_cuts_veto"        ,  "SF syst. unc. as Z p_{T} for Veto electron pairs (tag cuts)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_ele_generator_veto        = new TH1D("h_syst_Z_ele_generator_veto"       ,  "SF syst. unc. as Z p_{T} for Veto electron pairs (generator eff.)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_ele_background_shape_veto = new TH1D("h_syst_Z_ele_background_shape_veto",  "SF syst. unc. as Z p_{T} for Veto electron pairs (background shape)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_ele_signal_shape_veto     = new TH1D("h_syst_Z_ele_signal_shape_veto"    ,  "SF syst. unc. as Z p_{T} for Veto electron pairs (signal shape)",  nbins, 0, max_pT);
  
  TH1D *h_syst_Z_ele_combined_tight         = new TH1D("h_syst_Z_ele_combined_tight"        ,  "SF syst. unc. as Z p_{T} for Tight electron pairs (combined)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_ele_tag_cuts_tight         = new TH1D("h_syst_Z_ele_tag_cuts_tight"        ,  "SF syst. unc. as Z p_{T} for Tight electron pairs (tag cuts)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_ele_generator_tight        = new TH1D("h_syst_Z_ele_generator_tight"       ,  "SF syst. unc. as Z p_{T} for Tight electron pairs (generator eff.)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_ele_background_shape_tight = new TH1D("h_syst_Z_ele_background_shape_tight",  "SF syst. unc. as Z p_{T} for Tight electron pairs (background shape)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_ele_signal_shape_tight     = new TH1D("h_syst_Z_ele_signal_shape_tight"    ,  "SF syst. unc. as Z p_{T} for Tight electron pairs (signal shape)",  nbins, 0, max_pT);
  
  TH1D *h_syst_Z_mu_combined_loose         = new TH1D("h_syst_Z_mu_combined_loose"        ,  "SF syst. unc. as Z p_{T} for Loose muon pairs (combined)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_mu_tag_cuts_loose         = new TH1D("h_syst_Z_mu_tag_cuts_loose"        ,  "SF syst. unc. as Z p_{T} for Loose muon pairs (tag cuts)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_mu_generator_loose        = new TH1D("h_syst_Z_mu_generator_loose"       ,  "SF syst. unc. as Z p_{T} for Loose muon pairs (generator eff.)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_mu_background_shape_loose = new TH1D("h_syst_Z_mu_background_shape_loose",  "SF syst. unc. as Z p_{T} for Loose muon pairs (background shape)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_mu_signal_shape_loose     = new TH1D("h_syst_Z_mu_signal_shape_loose"    ,  "SF syst. unc. as Z p_{T} for Loose muon pairs (signal shape)",  nbins, 0, max_pT);
  
  TH1D *h_syst_Z_mu_combined_tight         = new TH1D("h_syst_Z_mu_combined_tight"        ,  "SF syst. unc. as Z p_{T} for Tight muon pairs (combined)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_mu_tag_cuts_tight         = new TH1D("h_syst_Z_mu_tag_cuts_tight"        ,  "SF syst. unc. as Z p_{T} for Tight muon pairs (tag cuts)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_mu_generator_tight        = new TH1D("h_syst_Z_mu_generator_tight"       ,  "SF syst. unc. as Z p_{T} for Tight muon pairs (generator eff.)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_mu_background_shape_tight = new TH1D("h_syst_Z_mu_background_shape_tight",  "SF syst. unc. as Z p_{T} for Tight muon pairs (background shape)",  nbins, 0, max_pT);
  TH1D *h_syst_Z_mu_signal_shape_tight     = new TH1D("h_syst_Z_mu_signal_shape_tight"    ,  "SF syst. unc. as Z p_{T} for Tight muon pairs (signal shape)",  nbins, 0, max_pT);
  
  TH1D *h_total_weight_ele_veto  = new TH1D("h_total_weight_ele_veto",  "",  nbins, 0, max_pT);
  TH1D *h_total_weight_ele_tight = new TH1D("h_total_weight_ele_tight", "", nbins, 0, max_pT);
  TH1D *h_total_weight_mu_loose  = new TH1D("h_total_weight_mu_loose",  "",  nbins, 0, max_pT);
  TH1D *h_total_weight_mu_tight  = new TH1D("h_total_weight_mu_tight",  "",  nbins, 0, max_pT);
  
  TTree *t_tnp_ele_veto = (TTree*) f_tnp_ele_veto->Get("Events");
  t_tnp_ele_veto->SetBranchAddress("runNum",   &runNum   );  
  t_tnp_ele_veto->SetBranchAddress("lumiSec",  &lumiSec  );  
  t_tnp_ele_veto->SetBranchAddress("evtNum",   &evtNum   );  
  t_tnp_ele_veto->SetBranchAddress("npv",      &npv      );  
  t_tnp_ele_veto->SetBranchAddress("pass",     &pass     );  
  t_tnp_ele_veto->SetBranchAddress("npu",      &npu      );  
  t_tnp_ele_veto->SetBranchAddress("scale1fb", &scale1fb );
  t_tnp_ele_veto->SetBranchAddress("mass",     &mass     );  
  t_tnp_ele_veto->SetBranchAddress("qtag",     &qtag     );  
  t_tnp_ele_veto->SetBranchAddress("qprobe",   &qprobe   );  
  t_tnp_ele_veto->SetBranchAddress("tag",      &p4_tag   );  
  t_tnp_ele_veto->SetBranchAddress("probe",    &p4_probe );     
  TTree *t_tnp_ele_tight = (TTree*) f_tnp_ele_tight->Get("Events");
  t_tnp_ele_tight->SetBranchAddress("runNum",   &runNum   );  
  t_tnp_ele_tight->SetBranchAddress("lumiSec",  &lumiSec  );  
  t_tnp_ele_tight->SetBranchAddress("evtNum",   &evtNum   );  
  t_tnp_ele_tight->SetBranchAddress("npv",      &npv      );  
  t_tnp_ele_tight->SetBranchAddress("pass",     &pass     );  
  t_tnp_ele_tight->SetBranchAddress("npu",      &npu      );  
  t_tnp_ele_tight->SetBranchAddress("scale1fb", &scale1fb );
  t_tnp_ele_tight->SetBranchAddress("mass",     &mass     );  
  t_tnp_ele_tight->SetBranchAddress("qtag",     &qtag     );  
  t_tnp_ele_tight->SetBranchAddress("qprobe",   &qprobe   );  
  t_tnp_ele_tight->SetBranchAddress("tag",      &p4_tag   );  
  t_tnp_ele_tight->SetBranchAddress("probe",    &p4_probe );     
  TTree *t_tnp_mu_loose = (TTree*) f_tnp_mu_loose->Get("Events");
  t_tnp_mu_loose->SetBranchAddress("runNum",   &runNum   );  
  t_tnp_mu_loose->SetBranchAddress("lumiSec",  &lumiSec  );  
  t_tnp_mu_loose->SetBranchAddress("evtNum",   &evtNum   );  
  t_tnp_mu_loose->SetBranchAddress("npv",      &npv      );  
  t_tnp_mu_loose->SetBranchAddress("pass",     &pass     );  
  t_tnp_mu_loose->SetBranchAddress("npu",      &npu      );  
  t_tnp_mu_loose->SetBranchAddress("scale1fb", &scale1fb );
  t_tnp_mu_loose->SetBranchAddress("mass",     &mass     );  
  t_tnp_mu_loose->SetBranchAddress("qtag",     &qtag     );  
  t_tnp_mu_loose->SetBranchAddress("qprobe",   &qprobe   );  
  t_tnp_mu_loose->SetBranchAddress("tag",      &p4_tag   );  
  t_tnp_mu_loose->SetBranchAddress("probe",    &p4_probe );     
  TTree *t_tnp_mu_tight = (TTree*) f_tnp_mu_tight->Get("Events");
  t_tnp_mu_tight->SetBranchAddress("runNum",   &runNum   );  
  t_tnp_mu_tight->SetBranchAddress("lumiSec",  &lumiSec  );  
  t_tnp_mu_tight->SetBranchAddress("evtNum",   &evtNum   );  
  t_tnp_mu_tight->SetBranchAddress("npv",      &npv      );  
  t_tnp_mu_tight->SetBranchAddress("pass",     &pass     );  
  t_tnp_mu_tight->SetBranchAddress("npu",      &npu      );  
  t_tnp_mu_tight->SetBranchAddress("scale1fb", &scale1fb );
  t_tnp_mu_tight->SetBranchAddress("mass",     &mass     );  
  t_tnp_mu_tight->SetBranchAddress("qtag",     &qtag     );  
  t_tnp_mu_tight->SetBranchAddress("qprobe",   &qprobe   );  
  t_tnp_mu_tight->SetBranchAddress("tag",      &p4_tag   );  
  t_tnp_mu_tight->SetBranchAddress("probe",    &p4_probe );     
  Long64_t nentries;
  
  printf("Doing veto electrons . . .\n");
  nentries= t_tnp_ele_veto->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    t_tnp_ele_veto->GetEntry(i);
    if(!(
      pass==1 &&
      p4_tag->Pt() >= 30 &&
      TMath::Abs(p4_tag->Eta()) <= 2.1 &&
      TMath::Abs(mass - 90) <= 30 &&
      qtag + qprobe == 0
    )) continue;
    int tag_bin=syst_ele_sf_combined_veto->FindBin(TMath::Abs(p4_tag->Eta()), p4_tag->Pt());
    int probe_bin=syst_ele_sf_combined_veto->FindBin(TMath::Abs(p4_probe->Eta()), p4_probe->Pt());
    double syst_tag   = syst_ele_sf_combined_veto->GetBinContent(tag_bin);
    double syst_probe = syst_ele_sf_combined_veto->GetBinContent(probe_bin);
    double syst_Z = syst_tag+syst_probe;
    double weighted_syst = scale1fb*(syst_Z);
    if(verbose) printf("tag pt = %f, eta = %f\nprobe pt = %f, eta = %f\n", p4_tag->Pt(), p4_tag->Eta(), p4_probe->Pt(), p4_probe->Eta());
    if(verbose) printf("syst error for entry %lld is %f + %f = %f\n", i, syst_tag, syst_probe, syst_Z);
    if(verbose) printf("weighted syst error is %f\n", weighted_syst);

    TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
    double Z_pT = systemP4.Pt();
    if(verbose) printf("Z pT is %f\n", Z_pT);
    h_Z_pT_ele_veto->Fill(Z_pT, scale1fb);
    h_total_weight_ele_veto->Fill(Z_pT,scale1fb);
    h_syst_Z_ele_combined_veto->Fill(Z_pT, weighted_syst);
    h_syst_Z_ele_tag_cuts_veto        ->Fill(Z_pT, scale1fb*(syst_ele_sf_tag_cuts_veto ->GetBinContent(tag_bin)        + syst_ele_sf_tag_cuts_veto->GetBinContent(probe_bin)));
    h_syst_Z_ele_generator_veto       ->Fill(Z_pT, scale1fb*(syst_ele_sf_generator_veto->GetBinContent(tag_bin)        + syst_ele_sf_generator_veto->GetBinContent(probe_bin)));
    h_syst_Z_ele_background_shape_veto->Fill(Z_pT, scale1fb*(syst_ele_sf_background_shape_veto->GetBinContent(tag_bin) + syst_ele_sf_background_shape_veto->GetBinContent(probe_bin)));
    h_syst_Z_ele_signal_shape_veto    ->Fill(Z_pT, scale1fb*(syst_ele_sf_signal_shape_veto->GetBinContent(tag_bin)     + syst_ele_sf_signal_shape_veto->GetBinContent(probe_bin)));
  }
  h_Z_pT_ele_veto    ->Scale(1./h_total_weight_ele_veto->Integral());
  h_syst_Z_ele_combined_veto         ->Divide(h_total_weight_ele_veto);
  h_syst_Z_ele_tag_cuts_veto         ->Divide(h_total_weight_ele_veto);
  h_syst_Z_ele_generator_veto        ->Divide(h_total_weight_ele_veto);
  h_syst_Z_ele_background_shape_veto ->Divide(h_total_weight_ele_veto);
  h_syst_Z_ele_signal_shape_veto     ->Divide(h_total_weight_ele_veto);
  
  
  printf("Doing tight electrons . . .\n");
  TObjArray *dilepton_pt_correlation_ele_tight_ = new TObjArray(nbins);
  for(int j=0; j<nbins; j++) {
    char name[128], title[128];
    sprintf(name, "dilepton_pt_correlation_ele_tight_%d", j);
    sprintf(title, "Correlation of tag and probe pT (tight electrons, Z p_{T} [%d, %d])", (int) j*max_pT/nbins, (j+1)*max_pT/nbins);
    TH2D *dilepton_pt_correlation = new TH2D(name, title, n_ele_pt_bins, ele_pt_bins, n_ele_pt_bins, ele_pt_bins);
    dilepton_pt_correlation->GetXaxis()->SetTitle("probe p_{T} [GeV]");
    dilepton_pt_correlation->GetYaxis()->SetTitle("tag p_{T} [GeV]");
    dilepton_pt_correlation->GetYaxis()->SetRangeUser(30, ele_pt_bins[n_ele_pt_bins]);
    dilepton_pt_correlation->SetMaximum(1.);
    dilepton_pt_correlation->SetMinimum(-0.001);
    dilepton_pt_correlation_ele_tight_->Add(dilepton_pt_correlation);
  }

  nentries= t_tnp_ele_tight->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    t_tnp_ele_tight->GetEntry(i);
    if(!(
      pass==1 &&
      p4_tag->Pt() >= 30 &&
      TMath::Abs(p4_tag->Eta()) <= 2.1 &&
      TMath::Abs(mass - 90) <= 30 &&
      qtag + qprobe == 0
    )) continue;
    int tag_bin=syst_ele_sf_combined_tight->FindBin(TMath::Abs(p4_tag->Eta()), p4_tag->Pt());
    int probe_bin=syst_ele_sf_combined_tight->FindBin(TMath::Abs(p4_probe->Eta()), p4_probe->Pt());
    double syst_tag   = syst_ele_sf_combined_tight->GetBinContent(syst_ele_sf_combined_tight->FindBin(TMath::Abs(p4_tag->Eta()), p4_tag->Pt()));
    double syst_probe = syst_ele_sf_combined_tight->GetBinContent(syst_ele_sf_combined_tight->FindBin(TMath::Abs(p4_probe->Eta()), p4_probe->Pt()));
    double syst_Z = syst_tag+syst_probe;
    double weighted_syst = scale1fb*(syst_Z);
    if(verbose) printf("tag pt = %f, eta = %f\nprobe pt = %f, eta = %f\n", p4_tag->Pt(), p4_tag->Eta(), p4_probe->Pt(), p4_probe->Eta());
    if(verbose) printf("syst error for entry %lld is %f + %f = %f\n", i, syst_tag, syst_probe, syst_Z);
    if(verbose) printf("weighted syst error is %f\n", weighted_syst);

    TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
    double Z_pT = systemP4.Pt();
    if(verbose) printf("Z pT is %f\n", Z_pT);
    h_Z_pT_ele_tight->Fill(Z_pT, scale1fb);
    h_total_weight_ele_tight->Fill(Z_pT,scale1fb);
    h_syst_Z_ele_combined_tight->Fill(Z_pT, weighted_syst);
    h_syst_Z_ele_tag_cuts_tight        ->Fill(Z_pT, scale1fb*(syst_ele_sf_tag_cuts_tight ->GetBinContent(tag_bin)        + syst_ele_sf_tag_cuts_tight->GetBinContent(probe_bin)));
    h_syst_Z_ele_generator_tight       ->Fill(Z_pT, scale1fb*(syst_ele_sf_generator_tight->GetBinContent(tag_bin)        + syst_ele_sf_generator_tight->GetBinContent(probe_bin)));
    h_syst_Z_ele_background_shape_tight->Fill(Z_pT, scale1fb*(syst_ele_sf_background_shape_tight->GetBinContent(tag_bin) + syst_ele_sf_background_shape_tight->GetBinContent(probe_bin)));
    h_syst_Z_ele_signal_shape_tight    ->Fill(Z_pT, scale1fb*(syst_ele_sf_signal_shape_tight->GetBinContent(tag_bin)     + syst_ele_sf_signal_shape_tight->GetBinContent(probe_bin)));
  }
  h_Z_pT_ele_tight    ->Scale(1./h_total_weight_ele_tight->Integral());
  h_syst_Z_ele_combined_tight         ->Divide(h_total_weight_ele_tight);
  h_syst_Z_ele_tag_cuts_tight         ->Divide(h_total_weight_ele_tight);
  h_syst_Z_ele_generator_tight        ->Divide(h_total_weight_ele_tight);
  h_syst_Z_ele_background_shape_tight ->Divide(h_total_weight_ele_tight);
  h_syst_Z_ele_signal_shape_tight     ->Divide(h_total_weight_ele_tight);

  printf("Doing loose muons . . .\n");
  nentries= t_tnp_mu_loose->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    t_tnp_mu_loose->GetEntry(i);
    if(!(
      pass==1 &&
      p4_tag->Pt() >= 30 &&
      TMath::Abs(p4_tag->Eta()) <= 2.1 &&
      TMath::Abs(mass - 90) <= 30 &&
      qtag + qprobe == 0
    )) continue;
    double tag_pt   = TMath::Min(119.9, p4_tag->Pt());
    double probe_pt = TMath::Min(119.9, p4_probe->Pt());
    double tag_abseta   = TMath::Abs(p4_tag->Eta());
    double probe_abseta = TMath::Abs(p4_probe->Eta());
    int tag_bin=syst_ele_sf_combined_tight->FindBin(TMath::Abs(p4_tag->Eta()), p4_tag->Pt());
    int probe_bin=syst_ele_sf_combined_tight->FindBin(TMath::Abs(p4_probe->Eta()), p4_probe->Pt());
    double syst_tag   = syst_mu_sf_combined_loose->GetBinContent(syst_mu_sf_combined_loose->FindBin(tag_abseta, tag_pt));
    double syst_probe = syst_mu_sf_combined_loose->GetBinContent(syst_mu_sf_combined_loose->FindBin(probe_abseta, probe_pt));
    double syst_Z = syst_tag+syst_probe;
    double weighted_syst = scale1fb*(syst_Z);
    if(verbose) printf("tag pt = %f, eta = %f\nprobe pt = %f, eta = %f\n", p4_tag->Pt(), p4_tag->Eta(), p4_probe->Pt(), p4_probe->Eta());
    if(verbose) printf("syst error for entry %lld is %f + %f = %f\n", i, syst_tag, syst_probe, syst_Z);
    if(verbose) printf("weighted syst error is %f\n", weighted_syst);

    TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
    double Z_pT = systemP4.Pt();
    if(verbose) printf("Z pT is %f\n", Z_pT);
    h_Z_pT_mu_loose->Fill(Z_pT, scale1fb);
    h_total_weight_mu_loose->Fill(Z_pT,scale1fb);
    h_syst_Z_mu_combined_loose->Fill(Z_pT, weighted_syst);
    h_syst_Z_mu_tag_cuts_loose        ->Fill(Z_pT, scale1fb*(syst_mu_sf_tag_cuts_loose ->GetBinContent(tag_bin)        + syst_mu_sf_tag_cuts_loose->GetBinContent(probe_bin)));
    h_syst_Z_mu_generator_loose       ->Fill(Z_pT, scale1fb*(syst_mu_sf_generator_loose->GetBinContent(tag_bin)        + syst_mu_sf_generator_loose->GetBinContent(probe_bin)));
    h_syst_Z_mu_background_shape_loose->Fill(Z_pT, scale1fb*(syst_mu_sf_background_shape_loose->GetBinContent(tag_bin) + syst_mu_sf_background_shape_loose->GetBinContent(probe_bin)));
    h_syst_Z_mu_signal_shape_loose    ->Fill(Z_pT, scale1fb*(syst_mu_sf_signal_shape_loose->GetBinContent(tag_bin)     + syst_mu_sf_signal_shape_loose->GetBinContent(probe_bin)));
  }
  h_Z_pT_mu_loose    ->Scale(1./h_total_weight_mu_loose->Integral());
  h_syst_Z_mu_combined_loose         ->Divide(h_total_weight_mu_loose);
  h_syst_Z_mu_tag_cuts_loose         ->Divide(h_total_weight_mu_loose);
  h_syst_Z_mu_generator_loose        ->Divide(h_total_weight_mu_loose);
  h_syst_Z_mu_background_shape_loose ->Divide(h_total_weight_mu_loose);
  h_syst_Z_mu_signal_shape_loose     ->Divide(h_total_weight_mu_loose);
  
  printf("Doing tight muons . . .\n");
  nentries= t_tnp_mu_tight->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    t_tnp_mu_tight->GetEntry(i);
    if(!(
      pass==1 &&
      p4_tag->Pt() >= 30 &&
      TMath::Abs(p4_tag->Eta()) <= 2.1 &&
      TMath::Abs(mass - 90) <= 30 &&
      qtag + qprobe == 0
    )) continue;
    double tag_pt   = TMath::Min(119.9, p4_tag->Pt());
    double probe_pt = TMath::Min(119.9, p4_probe->Pt());
    double tag_abseta   = TMath::Abs(p4_tag->Eta());
    double probe_abseta = TMath::Abs(p4_probe->Eta());
    int tag_bin=syst_ele_sf_combined_tight->FindBin(TMath::Abs(p4_tag->Eta()), p4_tag->Pt());
    int probe_bin=syst_ele_sf_combined_tight->FindBin(TMath::Abs(p4_probe->Eta()), p4_probe->Pt());
    double syst_tag   = syst_mu_sf_combined_tight->GetBinContent(syst_mu_sf_combined_tight->FindBin(tag_abseta, tag_pt));
    double syst_probe = syst_mu_sf_combined_tight->GetBinContent(syst_mu_sf_combined_tight->FindBin(probe_abseta, probe_pt));
    double syst_Z = syst_tag+syst_probe;
    double weighted_syst = scale1fb*(syst_Z);
   
    TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
    double Z_pT = systemP4.Pt();
    if(verbose) printf("Z pT is %f\n", Z_pT);
    h_Z_pT_mu_tight->Fill(Z_pT, scale1fb);
    h_total_weight_mu_tight->Fill(Z_pT,scale1fb);
    h_syst_Z_mu_combined_tight->Fill(Z_pT, weighted_syst);
    h_syst_Z_mu_tag_cuts_tight        ->Fill(Z_pT, scale1fb*(syst_mu_sf_tag_cuts_tight ->GetBinContent(tag_bin)        + syst_mu_sf_tag_cuts_tight->GetBinContent(probe_bin)));
    h_syst_Z_mu_generator_tight       ->Fill(Z_pT, scale1fb*(syst_mu_sf_generator_tight->GetBinContent(tag_bin)        + syst_mu_sf_generator_tight->GetBinContent(probe_bin)));
    h_syst_Z_mu_background_shape_tight->Fill(Z_pT, scale1fb*(syst_mu_sf_background_shape_tight->GetBinContent(tag_bin) + syst_mu_sf_background_shape_tight->GetBinContent(probe_bin)));
    h_syst_Z_mu_signal_shape_tight    ->Fill(Z_pT, scale1fb*(syst_mu_sf_signal_shape_tight->GetBinContent(tag_bin)     + syst_mu_sf_signal_shape_tight->GetBinContent(probe_bin)));
  }
  h_Z_pT_mu_tight    ->Scale(1./h_total_weight_mu_tight->Integral());
  h_syst_Z_mu_combined_tight         ->Divide(h_total_weight_mu_tight);
  h_syst_Z_mu_tag_cuts_tight         ->Divide(h_total_weight_mu_tight);
  h_syst_Z_mu_generator_tight        ->Divide(h_total_weight_mu_tight);
  h_syst_Z_mu_background_shape_tight ->Divide(h_total_weight_mu_tight);
  h_syst_Z_mu_signal_shape_tight     ->Divide(h_total_weight_mu_tight);
  
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);
  
  TBox *box_ele = new TBox(0,0,200,.1);
  box_ele->SetFillColor(mit_red);
  box_ele->SetFillStyle(3354);
  box_ele->SetLineColor(2);
  box_ele->SetLineStyle(1);
  box_ele->SetLineWidth(2);
  TBox *box_mu = new TBox(0,0,200,.1);
  box_mu->SetFillColor(mit_red);
  box_mu->SetFillStyle(3354);
  box_mu->SetLineColor(2);
  box_mu->SetLineStyle(1);
  box_mu->SetLineWidth(2);
  for(int j=1; j<=nbins; j++) {
    h_syst_Z_ele_combined_veto->SetBinError(j,0);
    h_syst_Z_ele_combined_tight->SetBinError(j,0);
    h_syst_Z_mu_combined_loose->SetBinError(j,0);
    h_syst_Z_mu_combined_tight->SetBinError(j,0);
  }

  TCanvas *c_Z_pT_ele_veto = new TCanvas("c_Z_pT_ele_veto","Z pT spectrum for Veto electron pairs");
  h_Z_pT_ele_veto->Draw();
  h_Z_pT_ele_veto->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  h_Z_pT_ele_veto->GetXaxis()->SetTitleOffset(1.25);
  gPad->Update();
  TCanvas *c_syst_Z_ele_combined_veto = new TCanvas("c_syst_Z_ele_combined_veto", "SF syst. unc. as Z pT for Veto electron pairs");
  c_syst_Z_ele_combined_veto->SetGrid();
  h_syst_Z_ele_combined_veto->SetMarkerStyle(21);
  h_syst_Z_ele_combined_veto->SetMarkerSize(1);
  h_syst_Z_ele_combined_veto->SetMarkerColor(1);
  h_syst_Z_ele_combined_veto->SetLineColor(1);
  h_syst_Z_ele_combined_veto->SetLineWidth(2);
  h_syst_Z_ele_combined_veto->SetMaximum(.1);
  h_syst_Z_ele_combined_veto->SetMinimum(0);
  h_syst_Z_ele_combined_veto->Draw("LP");
  h_syst_Z_ele_combined_veto->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  h_syst_Z_ele_combined_veto->GetYaxis()->SetTitle("Uncertainty");
  h_syst_Z_ele_combined_veto->GetXaxis()->SetTitleOffset(1.25);
  h_syst_Z_ele_combined_veto->GetYaxis()->SetTitleOffset(1.25);
  box_ele->Draw("L");
  h_syst_Z_ele_combined_veto->Draw("LP SAME");
  
  TCanvas *c_Z_pT_ele_tight = new TCanvas("c_Z_pT_ele_tight","Z pT spectrum for Tight electron pairs");
  h_Z_pT_ele_tight->Draw();
  h_Z_pT_ele_tight->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  h_Z_pT_ele_tight->GetXaxis()->SetTitleOffset(1.25);
  gPad->Update();
  TCanvas *c_syst_Z_ele_combined_tight = new TCanvas("c_syst_Z_ele_combined_tight", "SF syst. unc. as Z pT for Tight electron pairs");
  c_syst_Z_ele_combined_tight->SetGrid();
  h_syst_Z_ele_combined_tight->SetMarkerStyle(21);
  h_syst_Z_ele_combined_tight->SetMarkerSize(1);
  h_syst_Z_ele_combined_tight->SetMarkerColor(1);
  h_syst_Z_ele_combined_tight->SetLineColor(1);
  h_syst_Z_ele_combined_tight->SetLineWidth(2);
  h_syst_Z_ele_combined_tight->SetMaximum(.1);
  h_syst_Z_ele_combined_tight->SetMinimum(0);
  h_syst_Z_ele_combined_tight->Draw("LP");
  h_syst_Z_ele_combined_tight->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  h_syst_Z_ele_combined_tight->GetYaxis()->SetTitle("Uncertainty");
  h_syst_Z_ele_combined_tight->GetXaxis()->SetTitleOffset(1.25);
  h_syst_Z_ele_combined_tight->GetYaxis()->SetTitleOffset(1.25);
  box_ele->Draw("L");
  h_syst_Z_ele_combined_tight->Draw("LP SAME");

  
  TCanvas *c_Z_pT_mu_loose = new TCanvas("c_Z_pT_mu_loose","Z pT spectrum for Loose muon pairs");
  h_Z_pT_mu_loose->Draw();
  h_Z_pT_mu_loose->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  h_Z_pT_mu_loose->GetXaxis()->SetTitleOffset(1.25);
  gPad->Update();
  TCanvas *c_syst_Z_mu_combined_loose = new TCanvas("c_syst_Z_mu_combined_loose", "SF syst. unc. as Z pT for Loose muon pairs");
  c_syst_Z_mu_combined_loose->SetGrid();
  h_syst_Z_mu_combined_loose->SetMarkerStyle(21);
  h_syst_Z_mu_combined_loose->SetMarkerSize(1);
  h_syst_Z_mu_combined_loose->SetMarkerColor(1);
  h_syst_Z_mu_combined_loose->SetLineColor(1);
  h_syst_Z_mu_combined_loose->SetLineWidth(2);
  h_syst_Z_mu_combined_loose->SetMaximum(.1);
  h_syst_Z_mu_combined_loose->SetMinimum(0);
  h_syst_Z_mu_combined_loose->Draw("LP");
  h_syst_Z_mu_combined_loose->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  h_syst_Z_mu_combined_loose->GetYaxis()->SetTitle("Uncertainty");
  h_syst_Z_mu_combined_loose->GetXaxis()->SetTitleOffset(1.25);
  h_syst_Z_mu_combined_loose->GetYaxis()->SetTitleOffset(1.25);
  box_mu->Draw("L");
  h_syst_Z_mu_combined_loose->Draw("LP SAME");

  
  TCanvas *c_Z_pT_mu_tight = new TCanvas("c_Z_pT_mu_tight","Z pT spectrum for Tight muon pairs");
  h_Z_pT_mu_tight->Draw();
  h_Z_pT_mu_tight->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  h_Z_pT_mu_tight->GetXaxis()->SetTitleOffset(1.25);
  gPad->Update();
  TCanvas *c_syst_Z_mu_combined_tight = new TCanvas("c_syst_Z_mu_combined_tight", "SF syst. unc. as Z pT for Tight muon pairs");
  c_syst_Z_mu_combined_tight->SetGrid();
  h_syst_Z_mu_combined_tight->SetMarkerStyle(21);
  h_syst_Z_mu_combined_tight->SetMarkerSize(1);
  h_syst_Z_mu_combined_tight->SetMarkerColor(1);
  h_syst_Z_mu_combined_tight->SetLineColor(1);
  h_syst_Z_mu_combined_tight->SetLineWidth(2);
  h_syst_Z_mu_combined_tight->SetMaximum(.1);
  h_syst_Z_mu_combined_tight->SetMinimum(0);
  h_syst_Z_mu_combined_tight->Draw("LP");
  h_syst_Z_mu_combined_tight->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  h_syst_Z_mu_combined_tight->GetYaxis()->SetTitle("Uncertainty");
  h_syst_Z_mu_combined_tight->GetXaxis()->SetTitleOffset(1.25);
  h_syst_Z_mu_combined_tight->GetYaxis()->SetTitleOffset(1.25);
  box_mu->Draw("L");
  h_syst_Z_mu_combined_tight->Draw("LP SAME");

  c_Z_pT_ele_veto->Print((plots_dir+"Z_pT_ele_veto.png").c_str());
  c_syst_Z_ele_combined_veto->Print((plots_dir+"syst_Z_ele_veto.png").c_str());
  c_Z_pT_ele_tight->Print((plots_dir+"Z_pT_ele_tight.png").c_str());
  c_syst_Z_ele_combined_tight->Print((plots_dir+"syst_Z_ele_tight.png").c_str());
  c_Z_pT_mu_loose->Print((plots_dir+"Z_pT_mu_loose.png").c_str());
  c_syst_Z_mu_combined_loose->Print((plots_dir+"syst_Z_mu_loose.png").c_str());
  c_Z_pT_mu_tight->Print((plots_dir+"Z_pT_mu_tight.png").c_str());
  c_syst_Z_mu_combined_tight->Print((plots_dir+"syst_Z_mu_tight.png").c_str());

  mitPalette();
  gStyle->SetPaintTextFormat("4.3f");
  for(int j=0; j<nbins; j++) {
    TCanvas *c_dilepton_pt_correlation = new TCanvas;
    c_dilepton_pt_correlation->SetLogx();
    c_dilepton_pt_correlation->SetLogy();
    ((TH2D*) dilepton_pt_correlation_ele_tight_->At(j) )->Scale(1./((TH2D*) dilepton_pt_correlation_ele_tight_->At(j))->Integral() );
    ((TH2D*) dilepton_pt_correlation_ele_tight_->At(j) )->Draw("COLZ TEXT");
    c_dilepton_pt_correlation->Print((plots_dir + ((TH2D*) dilepton_pt_correlation_ele_tight_->At(j) )->GetName()+".png").c_str());
  }

  TCanvas *c_syst_breakdown_Z_ele_veto = new TCanvas("c_syst_breakdown_Z_ele_Veto", "SF syst. unc. by source as Z pT for Veto ele pairs");
  THStack *syst_breakdown_Z_ele_veto = new THStack("syst_breakdown_Z_ele_veto", "SF syst. unc. as Z p_{T} for Veto electron pairs (breakdown)");
  h_syst_Z_ele_signal_shape_veto     ->SetLineColor(mit_red);
  h_syst_Z_ele_background_shape_veto ->SetLineColor(mit_gray);
  h_syst_Z_ele_tag_cuts_veto         ->SetLineColor(38);
  h_syst_Z_ele_generator_veto        ->SetLineColor(53);
  h_syst_Z_ele_signal_shape_veto     ->SetFillColor(mit_red);
  h_syst_Z_ele_background_shape_veto ->SetFillColor(mit_gray);
  h_syst_Z_ele_tag_cuts_veto         ->SetFillColor(38);
  h_syst_Z_ele_generator_veto        ->SetFillColor(53);
  syst_breakdown_Z_ele_veto->Add(h_syst_Z_ele_signal_shape_veto, "B");
  syst_breakdown_Z_ele_veto->Add(h_syst_Z_ele_background_shape_veto, "B");
  syst_breakdown_Z_ele_veto->Add(h_syst_Z_ele_tag_cuts_veto, "B");
  syst_breakdown_Z_ele_veto->Add(h_syst_Z_ele_generator_veto, "B");
  TLegend *legend_syst_breakdown_Z_ele_veto = new TLegend(.6,.6,.8,.8);
  legend_syst_breakdown_Z_ele_veto->SetFillColor(0);
  legend_syst_breakdown_Z_ele_veto->AddEntry(h_syst_Z_ele_generator_veto, "MC Generator #varepsilon", "f");
  legend_syst_breakdown_Z_ele_veto->AddEntry(h_syst_Z_ele_tag_cuts_veto, "Tag cuts", "f");
  legend_syst_breakdown_Z_ele_veto->AddEntry(h_syst_Z_ele_background_shape_veto, "Background shape", "f");
  legend_syst_breakdown_Z_ele_veto->AddEntry(h_syst_Z_ele_signal_shape_veto, "Signal shape", "f");
  legend_syst_breakdown_Z_ele_veto->AddEntry(h_syst_Z_ele_combined_veto, "Combined", "lp");
  syst_breakdown_Z_ele_veto->SetMinimum(0);
  syst_breakdown_Z_ele_veto->SetMaximum(0.1);
  syst_breakdown_Z_ele_veto->Draw("HIST");
  h_syst_Z_ele_combined_veto->Draw("LP SAME");
  legend_syst_breakdown_Z_ele_veto->Draw("SAME");
  c_syst_breakdown_Z_ele_veto->Print((plots_dir+"syst_breakdown_Z_ele_veto.png").c_str());

  TCanvas *c_syst_breakdown_Z_ele_tight = new TCanvas("c_syst_breakdown_Z_ele_Tight", "SF syst. unc. by source as Z pT for Tight ele pairs");
  THStack *syst_breakdown_Z_ele_tight = new THStack("syst_breakdown_Z_ele_tight", "SF syst. unc. as Z p_{T} for Tight electron pairs (breakdown)");
  h_syst_Z_ele_signal_shape_tight     ->SetLineColor(mit_red);
  h_syst_Z_ele_background_shape_tight ->SetLineColor(mit_gray);
  h_syst_Z_ele_tag_cuts_tight         ->SetLineColor(38);
  h_syst_Z_ele_generator_tight        ->SetLineColor(53);
  h_syst_Z_ele_signal_shape_tight     ->SetFillColor(mit_red);
  h_syst_Z_ele_background_shape_tight ->SetFillColor(mit_gray);
  h_syst_Z_ele_tag_cuts_tight         ->SetFillColor(38);
  h_syst_Z_ele_generator_tight        ->SetFillColor(53);
  syst_breakdown_Z_ele_tight->Add(h_syst_Z_ele_signal_shape_tight, "B");
  syst_breakdown_Z_ele_tight->Add(h_syst_Z_ele_background_shape_tight, "B");
  syst_breakdown_Z_ele_tight->Add(h_syst_Z_ele_tag_cuts_tight, "B");
  syst_breakdown_Z_ele_tight->Add(h_syst_Z_ele_generator_tight, "B");
  TLegend *legend_syst_breakdown_Z_ele_tight = new TLegend(.6,.6,.8,.8);
  legend_syst_breakdown_Z_ele_tight->SetFillColor(0);
  legend_syst_breakdown_Z_ele_tight->AddEntry(h_syst_Z_ele_generator_tight, "MC Generator #varepsilon", "f");
  legend_syst_breakdown_Z_ele_tight->AddEntry(h_syst_Z_ele_tag_cuts_tight, "Tag cuts", "f");
  legend_syst_breakdown_Z_ele_tight->AddEntry(h_syst_Z_ele_background_shape_tight, "Background shape", "f");
  legend_syst_breakdown_Z_ele_tight->AddEntry(h_syst_Z_ele_signal_shape_tight, "Signal shape", "f");
  legend_syst_breakdown_Z_ele_tight->AddEntry(h_syst_Z_ele_combined_tight, "Combined", "lp");
  syst_breakdown_Z_ele_tight->SetMinimum(0);
  syst_breakdown_Z_ele_tight->SetMaximum(0.1);
  syst_breakdown_Z_ele_tight->Draw("HIST");
  h_syst_Z_ele_combined_tight->Draw("LP SAME");
  legend_syst_breakdown_Z_ele_tight->Draw("SAME");
  c_syst_breakdown_Z_ele_tight->Print((plots_dir+"syst_breakdown_Z_ele_tight.png").c_str());

  TCanvas *c_syst_breakdown_Z_mu_loose = new TCanvas("c_syst_breakdown_Z_mu_Loose", "SF syst. unc. by source as Z pT for Loose mu pairs");
  THStack *syst_breakdown_Z_mu_loose = new THStack("syst_breakdown_Z_mu_loose", "SF syst. unc. as Z p_{T} for Loose muon pairs (breakdown)");
  h_syst_Z_mu_signal_shape_loose     ->SetLineColor(mit_red);
  h_syst_Z_mu_background_shape_loose ->SetLineColor(mit_gray);
  h_syst_Z_mu_tag_cuts_loose         ->SetLineColor(38);
  h_syst_Z_mu_generator_loose        ->SetLineColor(53);
  h_syst_Z_mu_signal_shape_loose     ->SetFillColor(mit_red);
  h_syst_Z_mu_background_shape_loose ->SetFillColor(mit_gray);
  h_syst_Z_mu_tag_cuts_loose         ->SetFillColor(38);
  h_syst_Z_mu_generator_loose        ->SetFillColor(53);
  syst_breakdown_Z_mu_loose->Add(h_syst_Z_mu_signal_shape_loose, "B");
  syst_breakdown_Z_mu_loose->Add(h_syst_Z_mu_background_shape_loose, "B");
  syst_breakdown_Z_mu_loose->Add(h_syst_Z_mu_tag_cuts_loose, "B");
  syst_breakdown_Z_mu_loose->Add(h_syst_Z_mu_generator_loose, "B");
  TLegend *legend_syst_breakdown_Z_mu_loose = new TLegend(.6,.6,.8,.8);
  legend_syst_breakdown_Z_mu_loose->SetFillColor(0);
  legend_syst_breakdown_Z_mu_loose->AddEntry(h_syst_Z_mu_generator_loose, "MC Generator #varepsilon", "f");
  legend_syst_breakdown_Z_mu_loose->AddEntry(h_syst_Z_mu_tag_cuts_loose, "Tag cuts", "f");
  legend_syst_breakdown_Z_mu_loose->AddEntry(h_syst_Z_mu_background_shape_loose, "Background shape", "f");
  legend_syst_breakdown_Z_mu_loose->AddEntry(h_syst_Z_mu_signal_shape_loose, "Signal shape", "f");
  legend_syst_breakdown_Z_mu_loose->AddEntry(h_syst_Z_mu_combined_loose, "Combined", "lp");
  syst_breakdown_Z_mu_loose->SetMinimum(0);
  syst_breakdown_Z_mu_loose->SetMaximum(0.1);
  syst_breakdown_Z_mu_loose->Draw("HIST");
  h_syst_Z_mu_combined_loose->Draw("LP SAME"); 
  legend_syst_breakdown_Z_mu_loose->Draw("SAME");
  c_syst_breakdown_Z_mu_loose->Print((plots_dir+"syst_breakdown_Z_mu_loose.png").c_str());

  TCanvas *c_syst_breakdown_Z_mu_tight = new TCanvas("c_syst_breakdown_Z_mu_Tight", "SF syst. unc. by source as Z pT for Tight mu pairs");
  THStack *syst_breakdown_Z_mu_tight = new THStack("syst_breakdown_Z_mu_tight", "SF syst. unc. as Z p_{T} for Tight muon pairs (breakdown)");
  h_syst_Z_mu_signal_shape_tight     ->SetLineColor(mit_red);
  h_syst_Z_mu_background_shape_tight ->SetLineColor(mit_gray);
  h_syst_Z_mu_tag_cuts_tight         ->SetLineColor(38);
  h_syst_Z_mu_generator_tight        ->SetLineColor(53);
  h_syst_Z_mu_signal_shape_tight     ->SetFillColor(mit_red);
  h_syst_Z_mu_background_shape_tight ->SetFillColor(mit_gray);
  h_syst_Z_mu_tag_cuts_tight         ->SetFillColor(38);
  h_syst_Z_mu_generator_tight        ->SetFillColor(53);
  syst_breakdown_Z_mu_tight->Add(h_syst_Z_mu_signal_shape_tight, "B");
  syst_breakdown_Z_mu_tight->Add(h_syst_Z_mu_background_shape_tight, "B");
  syst_breakdown_Z_mu_tight->Add(h_syst_Z_mu_tag_cuts_tight, "B");
  syst_breakdown_Z_mu_tight->Add(h_syst_Z_mu_generator_tight, "B");
  TLegend *legend_syst_breakdown_Z_mu_tight = new TLegend(.6,.6,.8,.8);
  legend_syst_breakdown_Z_mu_tight->SetFillColor(0);
  legend_syst_breakdown_Z_mu_tight->AddEntry(h_syst_Z_mu_generator_tight, "MC Generator #varepsilon", "f");
  legend_syst_breakdown_Z_mu_tight->AddEntry(h_syst_Z_mu_tag_cuts_tight, "Tag cuts", "f");
  legend_syst_breakdown_Z_mu_tight->AddEntry(h_syst_Z_mu_background_shape_tight, "Background shape", "f");
  legend_syst_breakdown_Z_mu_tight->AddEntry(h_syst_Z_mu_signal_shape_tight, "Signal shape", "f");
  legend_syst_breakdown_Z_mu_tight->AddEntry(h_syst_Z_mu_combined_tight, "Combined", "lp");
  syst_breakdown_Z_mu_tight->SetMinimum(0);
  syst_breakdown_Z_mu_tight->SetMaximum(0.1);
  syst_breakdown_Z_mu_tight->Draw("HIST");
  h_syst_Z_mu_combined_tight->Draw("LP SAME");
  legend_syst_breakdown_Z_mu_tight->Draw("SAME");
  c_syst_breakdown_Z_mu_tight->Print((plots_dir+"syst_breakdown_Z_mu_tight.png").c_str());
  
  // Print out tables of systematics
  printf("########################################################\n");
  printf("# Veto Electron SF systematic                          #\n");
  printf("########################################################\n");
  printf("ZpT<=   Sig.    Bkg.    Tag     Gen.    Combined\n");
  for(int j=1; j<=nbins; j++) {
    printf("%8d  %5.4f  %5.4f  %5.4f  %5.4f  %5.4f\n", j*max_pT/nbins, h_syst_Z_ele_signal_shape_veto->GetBinContent(j), h_syst_Z_ele_background_shape_veto->GetBinContent(j), h_syst_Z_ele_tag_cuts_veto->GetBinContent(j), h_syst_Z_ele_generator_veto->GetBinContent(j), h_syst_Z_ele_combined_veto->GetBinContent(j));
  }
 
  printf("########################################################\n");
  printf("# Tight Electron SF systematic                         #\n");
  printf("########################################################\n\n");
  printf("ZpT<=   Sig.    Bkg.    Tag     Gen.    Combined\n");
  for(int j=1; j<=nbins; j++) {
    printf("%8d  %5.4f  %5.4f  %5.4f  %5.4f  %5.4f\n", j*max_pT/nbins, h_syst_Z_ele_signal_shape_tight->GetBinContent(j), h_syst_Z_ele_background_shape_tight->GetBinContent(j), h_syst_Z_ele_tag_cuts_tight->GetBinContent(j), h_syst_Z_ele_generator_tight->GetBinContent(j), h_syst_Z_ele_combined_tight->GetBinContent(j));
  }
 
  printf("########################################################\n");
  printf("# Loose Muon SF systematic                             #\n");
  printf("########################################################\n");
  printf("ZpT<=   Sig.    Bkg.    Tag     Gen.    Combined\n");
  for(int j=1; j<=nbins; j++) {
    printf("%8d  %5.4f  %5.4f  %5.4f  %5.4f  %5.4f\n", j*max_pT/nbins, h_syst_Z_mu_signal_shape_loose->GetBinContent(j), h_syst_Z_mu_background_shape_loose->GetBinContent(j), h_syst_Z_mu_tag_cuts_loose->GetBinContent(j), h_syst_Z_mu_generator_loose->GetBinContent(j), h_syst_Z_mu_combined_loose->GetBinContent(j));
  }
 
  printf("########################################################\n");
  printf("# Tight Muon SF systematic                             #\n");
  printf("########################################################\n");
  printf("ZpT<=   Sig.    Bkg.    Tag     Gen.    Combined\n");
  for(int j=1; j<=nbins; j++) {
    printf("%8d  %5.4f  %5.4f  %5.4f  %5.4f  %5.4f\n", j*max_pT/nbins, h_syst_Z_mu_signal_shape_tight->GetBinContent(j), h_syst_Z_mu_background_shape_tight->GetBinContent(j), h_syst_Z_mu_tag_cuts_tight->GetBinContent(j), h_syst_Z_mu_generator_tight->GetBinContent(j), h_syst_Z_mu_combined_tight->GetBinContent(j));
  }
 


}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void uncertainties_no_fact(
  string plots_dir,
  string root_dir,
  bool draw = true
) {
  gStyle->SetOptStat(0); 
  //open files
  TFile *f_mu_sf_BWCBPlusVoigt_erfcexp       = TFile::Open("~/leptonScaleFactors/root_local/16-02-2016/BWCBPlusVoigt_erfcexp/scalefactors_mu.root","READ");
  TFile *f_mu_sf_template_erfcexp = TFile::Open("~/leptonScaleFactors/root_local/16-02-2016/template_erfcexp/scalefactors_mu.root","UPDATE");
  TFile *f_mu_sf_template_exp     = TFile::Open("~/leptonScaleFactors/root_local/16-02-2016/template_exp/scalefactors_mu.root","READ");
  TFile *f_ele_sf_BWCBPlusVoigt_erfcexp       = TFile::Open("~/leptonScaleFactors/root_local/16-02-2016/BWCBPlusVoigt_erfcexp/scalefactors_ele.root","READ");
  TFile *f_ele_sf_template_erfcexp = TFile::Open("~/leptonScaleFactors/root_local/16-02-2016/template_erfcexp/scalefactors_ele.root","UPDATE");
  TFile *f_ele_sf_template_exp     = TFile::Open("~/leptonScaleFactors/root_local/16-02-2016/template_exp/scalefactors_ele.root","READ");
  
  // systematic uncertainty from signal model choice
  TH2D *syst_ele_sf_signal_veto   = (TH2D*) f_ele_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Veto_ele");
  TH2D *syst_ele_sf_signal_loose  = (TH2D*) f_ele_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Loose_ele");
  TH2D *syst_ele_sf_signal_medium = (TH2D*) f_ele_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Medium_ele");
  TH2D *syst_ele_sf_signal_tight  = (TH2D*) f_ele_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Tight_ele");
  TH2D *syst_mu_sf_signal_veto   = (TH2D*) f_mu_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Veto_mu");
  TH2D *syst_mu_sf_signal_loose  = (TH2D*) f_mu_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Loose_mu");
  TH2D *syst_mu_sf_signal_medium = (TH2D*) f_mu_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Medium_mu");
  TH2D *syst_mu_sf_signal_tight  = (TH2D*) f_mu_sf_BWCBPlusVoigt_erfcexp->Get("absDiff_unfactorized_scalefactors_Tight_mu");
  
  // systematic uncertainty from background model choice
  TH2D *syst_ele_sf_background_veto   = (TH2D*) f_ele_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Veto_ele");
  TH2D *syst_ele_sf_background_loose  = (TH2D*) f_ele_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Loose_ele");
  TH2D *syst_ele_sf_background_medium = (TH2D*) f_ele_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Medium_ele");
  TH2D *syst_ele_sf_background_tight  = (TH2D*) f_ele_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Tight_ele");
  TH2D *syst_mu_sf_background_veto   = (TH2D*) f_mu_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Veto_mu");
  TH2D *syst_mu_sf_background_loose  = (TH2D*) f_mu_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Loose_mu");
  TH2D *syst_mu_sf_background_medium = (TH2D*) f_mu_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Medium_mu");
  TH2D *syst_mu_sf_background_tight  = (TH2D*) f_mu_sf_template_exp->Get("absDiff_unfactorized_scalefactors_Tight_mu");

  // combined systematic uncertainty
  TH2D *syst_ele_sf_combined_veto   = new TH2D("syst_ele_sf_combined_veto", "Combined syst. unc. on electron veto scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *syst_ele_sf_combined_loose  = new TH2D("syst_ele_sf_combined_loose", "Combined syst. unc. on electron loose scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *syst_ele_sf_combined_medium = new TH2D("syst_ele_sf_combined_medium", "Combined syst. unc. on electron medium scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *syst_ele_sf_combined_tight  = new TH2D("syst_ele_sf_combined_tight", "Combined syst. unc. on electron tight scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *syst_mu_sf_combined_veto   = new TH2D("syst_mu_sf_combined_veto", "Combined syst. unc. on muon veto scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  TH2D *syst_mu_sf_combined_loose  = new TH2D("syst_mu_sf_combined_loose", "Combined syst. unc. on muon loose scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  TH2D *syst_mu_sf_combined_medium = new TH2D("syst_mu_sf_combined_medium", "Combined syst. unc. on muon medium scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  TH2D *syst_mu_sf_combined_tight  = new TH2D("syst_mu_sf_combined_tight", "Combined syst. unc. on muon tight scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  
  // statistical uncertainty
  TH2D *ele_sf_veto        = (TH2D*) f_ele_sf_template_erfcexp ->Get("unfactorized_scalefactors_Veto_ele;1"); 
  TH2D *ele_sf_loose       = (TH2D*) f_ele_sf_template_erfcexp ->Get("unfactorized_scalefactors_Loose_ele;1"); 
  TH2D *ele_sf_medium      = (TH2D*) f_ele_sf_template_erfcexp ->Get("unfactorized_scalefactors_Medium_ele;1"); 
  TH2D *ele_sf_tight       = (TH2D*) f_ele_sf_template_erfcexp ->Get("unfactorized_scalefactors_Tight_ele;1"); 
  TH2D *mu_sf_veto         = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Veto_mu;1");
  TH2D *mu_sf_loose        = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Loose_mu;1");
  TH2D *mu_sf_medium       = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Medium_mu;1");
  TH2D *mu_sf_tight        = (TH2D*) f_mu_sf_template_erfcexp  ->Get("unfactorized_scalefactors_Tight_mu;1");
  TH2D *stat_ele_sf_veto   = new TH2D("stat_ele_sf_veto",   "Stat. unc. on electron veto scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *stat_ele_sf_loose  = new TH2D("stat_ele_sf_loose",  "Stat. unc. on electron loose scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *stat_ele_sf_medium = new TH2D("stat_ele_sf_medium", "Stat. unc. on electron medium scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *stat_ele_sf_tight  = new TH2D("stat_ele_sf_tight",  "Stat. unc. on electron tight scalefactor",  n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  TH2D *stat_mu_sf_veto    = new TH2D("stat_mu_sf_veto",  "Stat. unc. on muon veto scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  TH2D *stat_mu_sf_loose   = new TH2D("stat_mu_sf_loose",  "Stat. unc. on muon loose scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  TH2D *stat_mu_sf_medium  = new TH2D("stat_mu_sf_medium",  "Stat. unc. on muon medium scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  TH2D *stat_mu_sf_tight   = new TH2D("stat_mu_sf_tight",  "Stat. unc. on muon tight scalefactor",  n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  
  for(int i_eta = 1; i_eta <= n_ele_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_ele_pt_bins; i_pt++) {
    int n = syst_ele_sf_combined_veto->GetBin(i_eta, i_pt);
    // assign combined electron syst. unc.
    syst_ele_sf_combined_veto->SetBinContent(n, sqrt(
      pow(syst_ele_sf_signal_veto     ->GetBinContent(n), 2) + 
      pow(syst_ele_sf_background_veto ->GetBinContent(n), 2) 
    )); 
    syst_ele_sf_combined_loose->SetBinContent(n, sqrt(
      pow(syst_ele_sf_signal_loose     ->GetBinContent(n), 2) + 
      pow(syst_ele_sf_background_loose ->GetBinContent(n), 2) 
    )); 
    syst_ele_sf_combined_medium->SetBinContent(n, sqrt(
      pow(syst_ele_sf_signal_medium     ->GetBinContent(n), 2) + 
      pow(syst_ele_sf_background_medium ->GetBinContent(n), 2) 
    )); 
    syst_ele_sf_combined_tight->SetBinContent(n, sqrt(
      pow(syst_ele_sf_signal_tight     ->GetBinContent(n), 2) + 
      pow(syst_ele_sf_background_tight ->GetBinContent(n), 2) 
    )); 
    // pull out electron stat. unc
    stat_ele_sf_veto->SetBinContent(n, ele_sf_veto->GetBinError(n));
    stat_ele_sf_loose->SetBinContent(n, ele_sf_loose->GetBinError(n));
    stat_ele_sf_medium->SetBinContent(n, ele_sf_medium->GetBinError(n));
    stat_ele_sf_tight->SetBinContent(n, ele_sf_tight->GetBinError(n));

    ele_sf_veto->SetBinError(n, sqrt(
      pow(stat_ele_sf_veto->GetBinContent(n), 2) +
      pow(syst_ele_sf_combined_veto->GetBinContent(n) , 2)
    ));
    ele_sf_loose->SetBinError(n, sqrt(
      pow(stat_ele_sf_loose->GetBinContent(n), 2) +
      pow(syst_ele_sf_combined_loose->GetBinContent(n) , 2)
    ));
    ele_sf_medium->SetBinError(n, sqrt(
      pow(stat_ele_sf_medium->GetBinContent(n), 2) +
      pow(syst_ele_sf_combined_medium->GetBinContent(n) , 2)
    ));
    ele_sf_tight->SetBinError(n, sqrt(
      pow(stat_ele_sf_tight->GetBinContent(n), 2) +
      pow(syst_ele_sf_combined_tight->GetBinContent(n) , 2)
    ));
  }}
  for(int i_eta = 1; i_eta <= n_mu_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_mu_pt_bins; i_pt++) {
    int n = syst_mu_sf_signal_veto->GetBin(i_eta, i_pt);
    // assign combined muon syst. unc.
    syst_mu_sf_combined_veto->SetBinContent(n, sqrt(
      pow(syst_mu_sf_signal_veto     ->GetBinContent(n), 2) +
      pow(syst_mu_sf_background_veto ->GetBinContent(n), 2)     
    )); 
    syst_mu_sf_combined_loose->SetBinContent(n, sqrt(
      pow(syst_mu_sf_signal_loose     ->GetBinContent(n), 2) +
      pow(syst_mu_sf_background_loose ->GetBinContent(n), 2)     
    )); 
    syst_mu_sf_combined_medium->SetBinContent(n, sqrt(
      pow(syst_mu_sf_signal_medium     ->GetBinContent(n), 2) +
      pow(syst_mu_sf_background_medium ->GetBinContent(n), 2)     
    )); 
    syst_mu_sf_combined_tight->SetBinContent(n, sqrt(
      pow(syst_mu_sf_signal_tight     ->GetBinContent(n), 2) + 
      pow(syst_mu_sf_background_tight ->GetBinContent(n), 2) 
    )); 
    // pull out electron stat. unc
    stat_mu_sf_veto   ->SetBinContent(n, mu_sf_veto->GetBinError(n));
    stat_mu_sf_loose  ->SetBinContent(n, mu_sf_loose->GetBinError(n));
    stat_mu_sf_medium ->SetBinContent(n, mu_sf_medium->GetBinError(n));
    stat_mu_sf_tight  ->SetBinContent(n, mu_sf_tight->GetBinError(n));

    mu_sf_veto->SetBinError(n, sqrt(
      pow(stat_mu_sf_veto->GetBinContent(n), 2) +
      pow(syst_mu_sf_combined_veto->GetBinContent(n) , 2)
    ));
    mu_sf_loose->SetBinError(n, sqrt(
      pow(stat_mu_sf_loose->GetBinContent(n), 2) +
      pow(syst_mu_sf_combined_loose->GetBinContent(n) , 2)
    ));
    mu_sf_medium->SetBinError(n, sqrt(
      pow(stat_mu_sf_medium->GetBinContent(n), 2) +
      pow(syst_mu_sf_combined_medium->GetBinContent(n) , 2)
    ));
    mu_sf_tight->SetBinError(n, sqrt(
      pow(stat_mu_sf_tight->GetBinContent(n), 2) +
      pow(syst_mu_sf_combined_tight->GetBinContent(n) , 2)
    ));

  }}
  f_ele_sf_template_erfcexp->cd();
  ele_sf_veto->Write("unfactorized_scalefactors_Veto_ele");
  ele_sf_loose->Write("unfactorized_scalefactors_Loose_ele");
  ele_sf_medium->Write("unfactorized_scalefactors_Medium_ele");
  ele_sf_tight->Write("unfactorized_scalefactors_Tight_ele");
  f_mu_sf_template_erfcexp->cd();
  mu_sf_veto->Write();
  mu_sf_loose->Write();
  mu_sf_medium->Write();
  mu_sf_tight->Write();
  if(draw) {
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.3f");
    mitPalette();
    TPaletteAxis *palette_axis;
    TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
    TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);
    TLine *bin1line=new TLine(20,0,20,0.1);
    TLine *bin2line=new TLine(30,0,30,0.1);
    TLine *bin3line=new TLine(40,0,40,0.1);
    TLine *bin4line=new TLine(50,0,50,0.1);
    TLine *bin5line=new TLine(70,0,70,0.1);
    bin1line->SetLineColor(1);
    bin2line->SetLineColor(1);
    bin3line->SetLineColor(1);
    bin4line->SetLineColor(1);
    bin5line->SetLineColor(1);
    bin1line->SetLineStyle(3);
    bin2line->SetLineStyle(3);
    bin3line->SetLineStyle(3);
    bin4line->SetLineStyle(3);
    bin5line->SetLineStyle(3);

    // ele veto plots
    TCanvas *c_syst_ele_sf_signal_veto = new TCanvas("c_syst_ele_sf_signal_veto", "Syst. unc. on ele. veto SF (signal)", 400, 600);
    syst_ele_sf_signal_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_signal_veto->SetTitle("Syst. unc. on ele. veto SF from signal choice");
    syst_ele_sf_signal_veto->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_signal_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_veto->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_veto->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_signal_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_veto->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_veto->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_veto->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_signal_veto->SetMinimum(0);
    syst_ele_sf_signal_veto->SetMaximum(0.2);
    syst_ele_sf_signal_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_signal_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_signal_veto->Update();
    c_syst_ele_sf_signal_veto->Print((plots_dir+"syst_ele_sf_signal_veto.png").c_str());
    
    TCanvas *c_syst_ele_sf_background_veto = new TCanvas("c_syst_ele_sf_background_veto", "Syst. unc. on ele. veto SF (b.g.)", 400, 600);
    syst_ele_sf_background_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_background_veto->SetTitle("Syst. unc. on ele. veto SF from background choice");
    syst_ele_sf_background_veto->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_background_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_veto->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_veto->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_background_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_veto->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_veto->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_veto->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_background_veto->SetMinimum(0);
    syst_ele_sf_background_veto->SetMaximum(0.2);
    syst_ele_sf_background_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_background_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_background_veto->Update();
    c_syst_ele_sf_background_veto->Print((plots_dir+"syst_ele_sf_background_veto.png").c_str());
    
    TCanvas *c_syst_ele_sf_combined_veto = new TCanvas("c_syst_ele_sf_combined_veto", "Syst. unc. on ele. veto SF (combined)", 400, 600);
    syst_ele_sf_combined_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_combined_veto->SetTitle("Syst. unc. on ele. veto SF (combined)");
    syst_ele_sf_combined_veto->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_combined_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_veto->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_veto->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_combined_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_veto->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_veto->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_veto->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_combined_veto->SetMinimum(0);
    syst_ele_sf_combined_veto->SetMaximum(0.2);
    syst_ele_sf_combined_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_combined_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_combined_veto->Update();
    c_syst_ele_sf_combined_veto->Print((plots_dir+"syst_ele_sf_combined_veto.png").c_str());
    
    TCanvas *c_syst_ele_sf_veto_1d = new TCanvas("c_syst_ele_sf_veto_1d","Systematic uncertainty on veto electron scale factors",800,600);
    c_syst_ele_sf_veto_1d->SetGrid(0,1);
    TH1D *syst_ele_sf_veto_eta1 = syst_ele_sf_combined_veto->ProjectionY("syst_ele_sf_veto_eta1", 1, 1);
    TH1D *syst_ele_sf_veto_eta2 = syst_ele_sf_combined_veto->ProjectionY("syst_ele_sf_veto_eta2", 2, 2);
    syst_ele_sf_veto_eta1->Draw("P0 L");
    syst_ele_sf_veto_eta1->SetTitle("Systematic uncertainty on veto electron scale factors");
    syst_ele_sf_veto_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_veto_eta1->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_veto_eta1->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_veto_eta1->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_veto_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
    syst_ele_sf_veto_eta1->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_veto_eta1->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_veto_eta1->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_veto_eta1->SetMinimum(0);
    syst_ele_sf_veto_eta1->SetMaximum(0.1);
    syst_ele_sf_veto_eta1->SetMarkerStyle(20);
    syst_ele_sf_veto_eta1->SetMarkerColor(1);
    syst_ele_sf_veto_eta1->SetLineColor(1);
    syst_ele_sf_veto_eta2->SetMarkerStyle(21);
    syst_ele_sf_veto_eta2->SetMarkerColor(1861);
    syst_ele_sf_veto_eta2->SetLineColor(1861);
    syst_ele_sf_veto_eta2->Draw("P0 L SAME");
    bin1line->Draw("SAME");
    bin2line->Draw("SAME");
    bin3line->Draw("SAME");
    bin4line->Draw("SAME");
    bin5line->Draw("SAME");
    TLegend *legend_veto_1d=new TLegend(.7,.7,.85,.85);
    legend_veto_1d->AddEntry(syst_ele_sf_veto_eta1,"| #eta | < 1.479", "lp");
    legend_veto_1d->AddEntry(syst_ele_sf_veto_eta2,"| #eta | > 1.479", "lp");
    legend_veto_1d->SetFillColor(0);
    legend_veto_1d->Draw("SAME");
    c_syst_ele_sf_veto_1d->Update();
    c_syst_ele_sf_veto_1d->Print((plots_dir+"syst_ele_sf_veto_1d.png").c_str());
  
    TCanvas *c_stat_ele_sf_veto = new TCanvas("c_stat_ele_sf_veto", "Stat. unc. on ele. veto SF", 400, 600);
    stat_ele_sf_veto->Draw("TEXT COLZ");
    gPad->Update();
    stat_ele_sf_veto->SetTitle("Stat. unc. on ele. veto SF");
    stat_ele_sf_veto->GetXaxis()->SetTitle("| #eta |");
    stat_ele_sf_veto->GetXaxis()->SetTitleOffset(0.9);
    stat_ele_sf_veto->GetXaxis()->SetTitleSize(0.04);
    stat_ele_sf_veto->GetXaxis()->SetLabelSize(0.02);
    stat_ele_sf_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_ele_sf_veto->GetYaxis()->SetTitleOffset(0.9);
    stat_ele_sf_veto->GetYaxis()->SetTitleSize(0.04);
    stat_ele_sf_veto->GetYaxis()->SetLabelSize(0.02);
    stat_ele_sf_veto->GetYaxis()->SetRangeUser(10,100);
    stat_ele_sf_veto->SetMinimum(0);
    stat_ele_sf_veto->SetMaximum(0.2);
    stat_ele_sf_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_ele_sf_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_ele_sf_veto->Update();
    c_stat_ele_sf_veto->Print((plots_dir+"stat_ele_sf_veto.png").c_str());

    // ele loose plots
    TCanvas *c_syst_ele_sf_signal_loose = new TCanvas("c_syst_ele_sf_signal_loose", "Syst. unc. on ele. loose SF (signal)", 400, 600);
    syst_ele_sf_signal_loose->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_signal_loose->SetTitle("Syst. unc. on ele. loose SF from signal choice");
    syst_ele_sf_signal_loose->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_signal_loose->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_loose->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_loose->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_signal_loose->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_loose->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_loose->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_loose->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_signal_loose->SetMinimum(0);
    syst_ele_sf_signal_loose->SetMaximum(0.2);
    syst_ele_sf_signal_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_signal_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_signal_loose->Update();
    c_syst_ele_sf_signal_loose->Print((plots_dir+"syst_ele_sf_signal_loose.png").c_str());
    
    TCanvas *c_syst_ele_sf_background_loose = new TCanvas("c_syst_ele_sf_background_loose", "Syst. unc. on ele. loose SF (b.g.)", 400, 600);
    syst_ele_sf_background_loose->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_background_loose->SetTitle("Syst. unc. on ele. loose SF from background choice");
    syst_ele_sf_background_loose->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_background_loose->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_loose->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_loose->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_background_loose->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_loose->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_loose->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_loose->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_background_loose->SetMinimum(0);
    syst_ele_sf_background_loose->SetMaximum(0.2);
    syst_ele_sf_background_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_background_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_background_loose->Update();
    c_syst_ele_sf_background_loose->Print((plots_dir+"syst_ele_sf_background_loose.png").c_str());
    
    TCanvas *c_syst_ele_sf_combined_loose = new TCanvas("c_syst_ele_sf_combined_loose", "Syst. unc. on ele. loose SF (combined)", 400, 600);
    syst_ele_sf_combined_loose->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_combined_loose->SetTitle("Syst. unc. on ele. loose SF (combined)");
    syst_ele_sf_combined_loose->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_combined_loose->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_loose->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_loose->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_combined_loose->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_loose->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_loose->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_loose->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_combined_loose->SetMinimum(0);
    syst_ele_sf_combined_loose->SetMaximum(0.2);
    syst_ele_sf_combined_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_combined_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_combined_loose->Update();
    c_syst_ele_sf_combined_loose->Print((plots_dir+"syst_ele_sf_combined_loose.png").c_str());
    
    TCanvas *c_syst_ele_sf_loose_1d = new TCanvas("c_syst_ele_sf_loose_1d","Systematic uncertainty on loose electron scale factors",800,600);
    c_syst_ele_sf_loose_1d->SetGrid(0,1);
    TH1D *syst_ele_sf_loose_eta1 = syst_ele_sf_combined_loose->ProjectionY("syst_ele_sf_loose_eta1", 1, 1);
    TH1D *syst_ele_sf_loose_eta2 = syst_ele_sf_combined_loose->ProjectionY("syst_ele_sf_loose_eta2", 2, 2);
    syst_ele_sf_loose_eta1->Draw("P0 L");
    syst_ele_sf_loose_eta1->SetTitle("Systematic uncertainty on loose electron scale factors");
    syst_ele_sf_loose_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_loose_eta1->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_loose_eta1->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_loose_eta1->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_loose_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
    syst_ele_sf_loose_eta1->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_loose_eta1->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_loose_eta1->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_loose_eta1->SetMinimum(0);
    syst_ele_sf_loose_eta1->SetMaximum(0.1);
    syst_ele_sf_loose_eta1->SetMarkerStyle(20);
    syst_ele_sf_loose_eta1->SetMarkerColor(1);
    syst_ele_sf_loose_eta1->SetLineColor(1);
    syst_ele_sf_loose_eta2->SetMarkerStyle(21);
    syst_ele_sf_loose_eta2->SetMarkerColor(1861);
    syst_ele_sf_loose_eta2->SetLineColor(1861);
    syst_ele_sf_loose_eta2->Draw("P0 L SAME");
    bin1line->Draw("SAME");
    bin2line->Draw("SAME");
    bin3line->Draw("SAME");
    bin4line->Draw("SAME");
    bin5line->Draw("SAME");
    TLegend *legend_loose_1d=new TLegend(.7,.7,.85,.85);
    legend_loose_1d->AddEntry(syst_ele_sf_loose_eta1,"| #eta | < 1.479", "lp");
    legend_loose_1d->AddEntry(syst_ele_sf_loose_eta2,"| #eta | > 1.479", "lp");
    legend_loose_1d->SetFillColor(0);
    legend_loose_1d->Draw("SAME");
    c_syst_ele_sf_loose_1d->Update();
    c_syst_ele_sf_loose_1d->Print((plots_dir+"syst_ele_sf_loose_1d.png").c_str());
   
    TCanvas *c_stat_ele_sf_loose = new TCanvas("c_stat_ele_sf_loose", "Stat. unc. on ele. loose SF", 400, 600);
    stat_ele_sf_loose->Draw("TEXT COLZ");
    gPad->Update();
    stat_ele_sf_loose->SetTitle("Stat. unc. on ele. loose SF");
    stat_ele_sf_loose->GetXaxis()->SetTitle("| #eta |");
    stat_ele_sf_loose->GetXaxis()->SetTitleOffset(0.9);
    stat_ele_sf_loose->GetXaxis()->SetTitleSize(0.04);
    stat_ele_sf_loose->GetXaxis()->SetLabelSize(0.02);
    stat_ele_sf_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_ele_sf_loose->GetYaxis()->SetTitleOffset(0.9);
    stat_ele_sf_loose->GetYaxis()->SetTitleSize(0.04);
    stat_ele_sf_loose->GetYaxis()->SetLabelSize(0.02);
    stat_ele_sf_loose->GetYaxis()->SetRangeUser(10,100);
    stat_ele_sf_loose->SetMinimum(0);
    stat_ele_sf_loose->SetMaximum(0.2);
    stat_ele_sf_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_ele_sf_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_ele_sf_loose->Update();
    c_stat_ele_sf_loose->Print((plots_dir+"stat_ele_sf_loose.png").c_str());

    // ele medium plots
    TCanvas *c_syst_ele_sf_signal_medium = new TCanvas("c_syst_ele_sf_signal_medium", "Syst. unc. on ele. medium SF (signal)", 400, 600);
    syst_ele_sf_signal_medium->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_signal_medium->SetTitle("Syst. unc. on ele. medium SF from signal choice");
    syst_ele_sf_signal_medium->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_signal_medium->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_medium->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_medium->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_medium->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_signal_medium->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_medium->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_medium->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_medium->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_signal_medium->SetMinimum(0);
    syst_ele_sf_signal_medium->SetMaximum(0.2);
    syst_ele_sf_signal_medium->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_signal_medium->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_signal_medium->Update();
    c_syst_ele_sf_signal_medium->Print((plots_dir+"syst_ele_sf_signal_medium.png").c_str());
    
    TCanvas *c_syst_ele_sf_background_medium = new TCanvas("c_syst_ele_sf_background_medium", "Syst. unc. on ele. medium SF (b.g.)", 400, 600);
    syst_ele_sf_background_medium->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_background_medium->SetTitle("Syst. unc. on ele. medium SF from background choice");
    syst_ele_sf_background_medium->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_background_medium->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_medium->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_medium->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_medium->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_background_medium->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_medium->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_medium->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_medium->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_background_medium->SetMinimum(0);
    syst_ele_sf_background_medium->SetMaximum(0.2);
    syst_ele_sf_background_medium->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_background_medium->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_background_medium->Update();
    c_syst_ele_sf_background_medium->Print((plots_dir+"syst_ele_sf_background_medium.png").c_str());
    
    TCanvas *c_syst_ele_sf_combined_medium = new TCanvas("c_syst_ele_sf_combined_medium", "Syst. unc. on ele. medium SF (combined)", 400, 600);
    syst_ele_sf_combined_medium->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_combined_medium->SetTitle("Syst. unc. on ele. medium SF (combined)");
    syst_ele_sf_combined_medium->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_combined_medium->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_medium->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_medium->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_medium->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_combined_medium->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_medium->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_medium->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_medium->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_combined_medium->SetMinimum(0);
    syst_ele_sf_combined_medium->SetMaximum(0.2);
    syst_ele_sf_combined_medium->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_combined_medium->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_combined_medium->Update();
    c_syst_ele_sf_combined_medium->Print((plots_dir+"syst_ele_sf_combined_medium.png").c_str());
    
    TCanvas *c_syst_ele_sf_medium_1d = new TCanvas("c_syst_ele_sf_medium_1d","Systematic uncertainty on medium electron scale factors",800,600);
    c_syst_ele_sf_medium_1d->SetGrid(0,1);
    TH1D *syst_ele_sf_medium_eta1 = syst_ele_sf_combined_medium->ProjectionY("syst_ele_sf_medium_eta1", 1, 1);
    TH1D *syst_ele_sf_medium_eta2 = syst_ele_sf_combined_medium->ProjectionY("syst_ele_sf_medium_eta2", 2, 2);
    syst_ele_sf_medium_eta1->Draw("P0 L");
    syst_ele_sf_medium_eta1->SetTitle("Systematic uncertainty on medium electron scale factors");
    syst_ele_sf_medium_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_medium_eta1->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_medium_eta1->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_medium_eta1->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_medium_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
    syst_ele_sf_medium_eta1->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_medium_eta1->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_medium_eta1->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_medium_eta1->SetMinimum(0);
    syst_ele_sf_medium_eta1->SetMaximum(0.1);
    syst_ele_sf_medium_eta1->SetMarkerStyle(20);
    syst_ele_sf_medium_eta1->SetMarkerColor(1);
    syst_ele_sf_medium_eta1->SetLineColor(1);
    syst_ele_sf_medium_eta2->SetMarkerStyle(21);
    syst_ele_sf_medium_eta2->SetMarkerColor(1861);
    syst_ele_sf_medium_eta2->SetLineColor(1861);
    syst_ele_sf_medium_eta2->Draw("P0 L SAME");
    bin1line->Draw("SAME");
    bin2line->Draw("SAME");
    bin3line->Draw("SAME");
    bin4line->Draw("SAME");
    bin5line->Draw("SAME");
    TLegend *legend_medium_1d=new TLegend(.7,.7,.85,.85);
    legend_medium_1d->AddEntry(syst_ele_sf_medium_eta1,"| #eta | < 1.479", "lp");
    legend_medium_1d->AddEntry(syst_ele_sf_medium_eta2,"| #eta | > 1.479", "lp");
    legend_medium_1d->SetFillColor(0);
    legend_medium_1d->Draw("SAME");
    c_syst_ele_sf_medium_1d->Update();
    c_syst_ele_sf_medium_1d->Print((plots_dir+"syst_ele_sf_medium_1d.png").c_str());
  
    TCanvas *c_stat_ele_sf_medium = new TCanvas("c_stat_ele_sf_medium", "Stat. unc. on ele. medium SF", 400, 600);
    stat_ele_sf_medium->Draw("TEXT COLZ");
    gPad->Update();
    stat_ele_sf_medium->SetTitle("Stat. unc. on ele. medium SF");
    stat_ele_sf_medium->GetXaxis()->SetTitle("| #eta |");
    stat_ele_sf_medium->GetXaxis()->SetTitleOffset(0.9);
    stat_ele_sf_medium->GetXaxis()->SetTitleSize(0.04);
    stat_ele_sf_medium->GetXaxis()->SetLabelSize(0.02);
    stat_ele_sf_medium->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_ele_sf_medium->GetYaxis()->SetTitleOffset(0.9);
    stat_ele_sf_medium->GetYaxis()->SetTitleSize(0.04);
    stat_ele_sf_medium->GetYaxis()->SetLabelSize(0.02);
    stat_ele_sf_medium->GetYaxis()->SetRangeUser(10,100);
    stat_ele_sf_medium->SetMinimum(0);
    stat_ele_sf_medium->SetMaximum(0.2);
    stat_ele_sf_medium->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_ele_sf_medium->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_ele_sf_medium->Update();
    c_stat_ele_sf_medium->Print((plots_dir+"stat_ele_sf_medium.png").c_str());

    // ele tight plots
    TCanvas *c_syst_ele_sf_signal_tight = new TCanvas("c_syst_ele_sf_signal_tight", "Syst. unc. on ele. tight SF (signal)", 400, 600);
    syst_ele_sf_signal_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_signal_tight->SetTitle("Syst. unc. on ele. tight SF from signal choice");
    syst_ele_sf_signal_tight->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_signal_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_tight->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_tight->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_signal_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_signal_tight->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_signal_tight->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_signal_tight->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_signal_tight->SetMinimum(0);
    syst_ele_sf_signal_tight->SetMaximum(0.2);
    syst_ele_sf_signal_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_signal_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_signal_tight->Update();
    c_syst_ele_sf_signal_tight->Print((plots_dir+"syst_ele_sf_signal_tight.png").c_str());
    
    TCanvas *c_syst_ele_sf_background_tight = new TCanvas("c_syst_ele_sf_background_tight", "Syst. unc. on ele. tight SF (b.g.)", 400, 600);
    syst_ele_sf_background_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_background_tight->SetTitle("Syst. unc. on ele. tight SF from background choice");
    syst_ele_sf_background_tight->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_background_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_tight->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_tight->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_background_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_background_tight->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_background_tight->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_background_tight->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_background_tight->SetMinimum(0);
    syst_ele_sf_background_tight->SetMaximum(0.2);
    syst_ele_sf_background_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_background_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_background_tight->Update();
    c_syst_ele_sf_background_tight->Print((plots_dir+"syst_ele_sf_background_tight.png").c_str());
    
    TCanvas *c_syst_ele_sf_combined_tight = new TCanvas("c_syst_ele_sf_combined_tight", "Syst. unc. on ele. tight SF (combined)", 400, 600);
    syst_ele_sf_combined_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_ele_sf_combined_tight->SetTitle("Syst. unc. on ele. tight SF (combined)");
    syst_ele_sf_combined_tight->GetXaxis()->SetTitle("| #eta |");
    syst_ele_sf_combined_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_tight->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_tight->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_combined_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_combined_tight->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_combined_tight->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_combined_tight->GetYaxis()->SetRangeUser(10,100);
    syst_ele_sf_combined_tight->SetMinimum(0);
    syst_ele_sf_combined_tight->SetMaximum(0.2);
    syst_ele_sf_combined_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_ele_sf_combined_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_ele_sf_combined_tight->Update();
    c_syst_ele_sf_combined_tight->Print((plots_dir+"syst_ele_sf_combined_tight.png").c_str());
    
    TCanvas *c_syst_ele_sf_tight_1d = new TCanvas("c_syst_ele_sf_tight_1d","Systematic uncertainty on tight electron scale factors",800,600);
    c_syst_ele_sf_tight_1d->SetGrid(0,1);
    TH1D *syst_ele_sf_tight_eta1 = syst_ele_sf_combined_tight->ProjectionY("syst_ele_sf_tight_eta1", 1, 1);
    TH1D *syst_ele_sf_tight_eta2 = syst_ele_sf_combined_tight->ProjectionY("syst_ele_sf_tight_eta2", 2, 2);
    syst_ele_sf_tight_eta1->Draw("P0 L");
    syst_ele_sf_tight_eta1->SetTitle("Systematic uncertainty on tight electron scale factors");
    syst_ele_sf_tight_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
    syst_ele_sf_tight_eta1->GetXaxis()->SetTitleOffset(0.9);
    syst_ele_sf_tight_eta1->GetXaxis()->SetTitleSize(0.04);
    syst_ele_sf_tight_eta1->GetXaxis()->SetLabelSize(0.02);
    syst_ele_sf_tight_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
    syst_ele_sf_tight_eta1->GetYaxis()->SetTitleOffset(0.9);
    syst_ele_sf_tight_eta1->GetYaxis()->SetTitleSize(0.04);
    syst_ele_sf_tight_eta1->GetYaxis()->SetLabelSize(0.02);
    syst_ele_sf_tight_eta1->SetMinimum(0);
    syst_ele_sf_tight_eta1->SetMaximum(0.1);
    syst_ele_sf_tight_eta1->SetMarkerStyle(20);
    syst_ele_sf_tight_eta1->SetMarkerColor(1);
    syst_ele_sf_tight_eta1->SetLineColor(1);
    syst_ele_sf_tight_eta2->SetMarkerStyle(21);
    syst_ele_sf_tight_eta2->SetMarkerColor(1861);
    syst_ele_sf_tight_eta2->SetLineColor(1861);
    syst_ele_sf_tight_eta2->Draw("P0 L SAME");
    bin1line->Draw("SAME");
    bin2line->Draw("SAME");
    bin3line->Draw("SAME");
    bin4line->Draw("SAME");
    bin5line->Draw("SAME");
    TLegend *legend_tight_1d=new TLegend(.7,.7,.85,.85);
    legend_tight_1d->AddEntry(syst_ele_sf_tight_eta1,"| #eta | < 1.479", "lp");
    legend_tight_1d->AddEntry(syst_ele_sf_tight_eta2,"| #eta | > 1.479", "lp");
    legend_tight_1d->SetFillColor(0);
    legend_tight_1d->Draw("SAME");
    c_syst_ele_sf_tight_1d->Update();
    c_syst_ele_sf_tight_1d->Print((plots_dir+"syst_ele_sf_tight_1d.png").c_str());
 
    TCanvas *c_stat_ele_sf_tight = new TCanvas("c_stat_ele_sf_tight", "Stat. unc. on ele. tight SF", 400, 600);
    stat_ele_sf_tight->Draw("TEXT COLZ");
    gPad->Update();
    stat_ele_sf_tight->SetTitle("Stat. unc. on ele. tight SF");
    stat_ele_sf_tight->GetXaxis()->SetTitle("| #eta |");
    stat_ele_sf_tight->GetXaxis()->SetTitleOffset(0.9);
    stat_ele_sf_tight->GetXaxis()->SetTitleSize(0.04);
    stat_ele_sf_tight->GetXaxis()->SetLabelSize(0.02);
    stat_ele_sf_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_ele_sf_tight->GetYaxis()->SetTitleOffset(0.9);
    stat_ele_sf_tight->GetYaxis()->SetTitleSize(0.04);
    stat_ele_sf_tight->GetYaxis()->SetLabelSize(0.02);
    stat_ele_sf_tight->GetYaxis()->SetRangeUser(10,100);
    stat_ele_sf_tight->SetMinimum(0);
    stat_ele_sf_tight->SetMaximum(0.2);
    stat_ele_sf_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_ele_sf_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_ele_sf_tight->Update();
    c_stat_ele_sf_tight->Print((plots_dir+"stat_ele_sf_tight.png").c_str());

    // mu veto plots
    TCanvas *c_syst_mu_sf_signal_veto = new TCanvas("c_syst_mu_sf_signal_veto", "Syst. unc. on mu. veto SF (signal)", 400, 600);
    syst_mu_sf_signal_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_signal_veto->SetTitle("Syst. unc. on mu. veto SF from signal choice");
    syst_mu_sf_signal_veto->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_signal_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_veto->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_veto->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_signal_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_veto->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_veto->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_veto->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_signal_veto->SetMinimum(0);
    syst_mu_sf_signal_veto->SetMaximum(0.2);
    syst_mu_sf_signal_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_signal_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_signal_veto->Update();
    c_syst_mu_sf_signal_veto->Print((plots_dir+"syst_mu_sf_signal_veto.png").c_str());
    
    TCanvas *c_syst_mu_sf_background_veto = new TCanvas("c_syst_mu_sf_background_veto", "Syst. unc. on mu. veto SF (b.g.)", 400, 600);
    syst_mu_sf_background_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_background_veto->SetTitle("Syst. unc. on mu. veto SF from background choice");
    syst_mu_sf_background_veto->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_background_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_veto->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_veto->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_background_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_veto->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_veto->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_veto->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_background_veto->SetMinimum(0);
    syst_mu_sf_background_veto->SetMaximum(0.2);
    syst_mu_sf_background_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_background_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_background_veto->Update();
    c_syst_mu_sf_background_veto->Print((plots_dir+"syst_mu_sf_background_veto.png").c_str());
    
    TCanvas *c_syst_mu_sf_combined_veto = new TCanvas("c_syst_mu_sf_combined_veto", "Syst. unc. on mu. veto SF (combined)", 400, 600);
    syst_mu_sf_combined_veto->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_combined_veto->SetTitle("Syst. unc. on mu. veto SF (combined)");
    syst_mu_sf_combined_veto->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_combined_veto->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_veto->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_veto->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_combined_veto->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_veto->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_veto->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_veto->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_combined_veto->SetMinimum(0);
    syst_mu_sf_combined_veto->SetMaximum(0.2);
    syst_mu_sf_combined_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_combined_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_combined_veto->Update();
    c_syst_mu_sf_combined_veto->Print((plots_dir+"syst_mu_sf_combined_veto.png").c_str());
    
    TCanvas *c_syst_mu_sf_veto_1d = new TCanvas("c_syst_mu_sf_veto_1d","Systematic uncertainty on veto muon scale factors",800,600);
    c_syst_mu_sf_veto_1d->SetGrid(0,1);
    TH1D *syst_mu_sf_veto_eta1 = syst_mu_sf_combined_veto->ProjectionY("syst_mu_sf_veto_eta1", 1, 1);
    TH1D *syst_mu_sf_veto_eta2 = syst_mu_sf_combined_veto->ProjectionY("syst_mu_sf_veto_eta2", 2, 2);
    syst_mu_sf_veto_eta1->Draw("P0 L");
    syst_mu_sf_veto_eta1->SetTitle("Systematic uncertainty on veto muon scale factors");
    syst_mu_sf_veto_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_veto_eta1->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_veto_eta1->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_veto_eta1->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_veto_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
    syst_mu_sf_veto_eta1->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_veto_eta1->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_veto_eta1->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_veto_eta1->SetMinimum(0);
    syst_mu_sf_veto_eta1->SetMaximum(0.1);
    syst_mu_sf_veto_eta1->SetMarkerStyle(20);
    syst_mu_sf_veto_eta1->SetMarkerColor(1);
    syst_mu_sf_veto_eta1->SetLineColor(1);
    syst_mu_sf_veto_eta2->SetMarkerStyle(21);
    syst_mu_sf_veto_eta2->SetMarkerColor(1861);
    syst_mu_sf_veto_eta2->SetLineColor(1861);
    syst_mu_sf_veto_eta2->Draw("P0 L SAME");
    bin1line->Draw("SAME");
    bin2line->Draw("SAME");
    bin3line->Draw("SAME");
    bin4line->Draw("SAME");
    bin5line->Draw("SAME");
    TLegend *legend_mu_veto_1d=new TLegend(.7,.7,.85,.85);
    legend_mu_veto_1d->AddEntry(syst_mu_sf_veto_eta1,"| #eta | < 1.479", "lp");
    legend_mu_veto_1d->AddEntry(syst_mu_sf_veto_eta2,"| #eta | > 1.479", "lp");
    legend_mu_veto_1d->SetFillColor(0);
    legend_mu_veto_1d->Draw("SAME");
    c_syst_mu_sf_veto_1d->Update();
    c_syst_mu_sf_veto_1d->Print((plots_dir+"syst_mu_sf_veto_1d.png").c_str());

    TCanvas *c_stat_mu_sf_veto = new TCanvas("c_stat_mu_sf_veto", "Stat. unc. on mu. veto SF", 400, 600);
    stat_mu_sf_veto->Draw("TEXT COLZ");
    gPad->Update();
    stat_mu_sf_veto->SetTitle("Stat. unc. on mu. veto SF");
    stat_mu_sf_veto->GetXaxis()->SetTitle("| #eta |");
    stat_mu_sf_veto->GetXaxis()->SetTitleOffset(0.9);
    stat_mu_sf_veto->GetXaxis()->SetTitleSize(0.04);
    stat_mu_sf_veto->GetXaxis()->SetLabelSize(0.02);
    stat_mu_sf_veto->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_mu_sf_veto->GetYaxis()->SetTitleOffset(0.9);
    stat_mu_sf_veto->GetYaxis()->SetTitleSize(0.04);
    stat_mu_sf_veto->GetYaxis()->SetLabelSize(0.02);
    stat_mu_sf_veto->GetYaxis()->SetRangeUser(10,100);
    stat_mu_sf_veto->SetMinimum(0);
    stat_mu_sf_veto->SetMaximum(0.2);
    stat_mu_sf_veto->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_mu_sf_veto->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_mu_sf_veto->Update();
    c_stat_mu_sf_veto->Print((plots_dir+"stat_mu_sf_veto.png").c_str());

    // mu loose plots
    TCanvas *c_syst_mu_sf_signal_loose = new TCanvas("c_syst_mu_sf_signal_loose", "Syst. unc. on mu. loose SF (signal)", 400, 600);
    syst_mu_sf_signal_loose->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_signal_loose->SetTitle("Syst. unc. on mu. loose SF from signal choice");
    syst_mu_sf_signal_loose->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_signal_loose->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_loose->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_loose->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_signal_loose->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_loose->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_loose->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_loose->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_signal_loose->SetMinimum(0);
    syst_mu_sf_signal_loose->SetMaximum(0.2);
    syst_mu_sf_signal_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_signal_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_signal_loose->Update();
    c_syst_mu_sf_signal_loose->Print((plots_dir+"syst_mu_sf_signal_loose.png").c_str());
    
    TCanvas *c_syst_mu_sf_background_loose = new TCanvas("c_syst_mu_sf_background_loose", "Syst. unc. on mu. loose SF (b.g.)", 400, 600);
    syst_mu_sf_background_loose->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_background_loose->SetTitle("Syst. unc. on mu. loose SF from background choice");
    syst_mu_sf_background_loose->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_background_loose->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_loose->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_loose->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_background_loose->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_loose->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_loose->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_loose->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_background_loose->SetMinimum(0);
    syst_mu_sf_background_loose->SetMaximum(0.2);
    syst_mu_sf_background_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_background_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_background_loose->Update();
    c_syst_mu_sf_background_loose->Print((plots_dir+"syst_mu_sf_background_loose.png").c_str());
    
    TCanvas *c_syst_mu_sf_combined_loose = new TCanvas("c_syst_mu_sf_combined_loose", "Syst. unc. on mu. loose SF (combined)", 400, 600);
    syst_mu_sf_combined_loose->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_combined_loose->SetTitle("Syst. unc. on mu. loose SF (combined)");
    syst_mu_sf_combined_loose->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_combined_loose->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_loose->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_loose->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_combined_loose->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_loose->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_loose->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_loose->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_combined_loose->SetMinimum(0);
    syst_mu_sf_combined_loose->SetMaximum(0.2);
    syst_mu_sf_combined_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_combined_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_combined_loose->Update();
    c_syst_mu_sf_combined_loose->Print((plots_dir+"syst_mu_sf_combined_loose.png").c_str());
    
    TCanvas *c_syst_mu_sf_loose_1d = new TCanvas("c_syst_mu_sf_loose_1d","Systematic uncertainty on loose muon scale factors",800,600);
    c_syst_mu_sf_loose_1d->SetGrid(0,1);
    TH1D *syst_mu_sf_loose_eta1 = syst_mu_sf_combined_loose->ProjectionY("syst_mu_sf_loose_eta1", 1, 1);
    TH1D *syst_mu_sf_loose_eta2 = syst_mu_sf_combined_loose->ProjectionY("syst_mu_sf_loose_eta2", 2, 2);
    syst_mu_sf_loose_eta1->Draw("P0 L");
    syst_mu_sf_loose_eta1->SetTitle("Systematic uncertainty on loose muon scale factors");
    syst_mu_sf_loose_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_loose_eta1->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_loose_eta1->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_loose_eta1->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_loose_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
    syst_mu_sf_loose_eta1->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_loose_eta1->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_loose_eta1->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_loose_eta1->SetMinimum(0);
    syst_mu_sf_loose_eta1->SetMaximum(0.1);
    syst_mu_sf_loose_eta1->SetMarkerStyle(20);
    syst_mu_sf_loose_eta1->SetMarkerColor(1);
    syst_mu_sf_loose_eta1->SetLineColor(1);
    syst_mu_sf_loose_eta2->SetMarkerStyle(21);
    syst_mu_sf_loose_eta2->SetMarkerColor(1861);
    syst_mu_sf_loose_eta2->SetLineColor(1861);
    syst_mu_sf_loose_eta2->Draw("P0 L SAME");
    bin1line->Draw("SAME");
    bin2line->Draw("SAME");
    bin3line->Draw("SAME");
    bin4line->Draw("SAME");
    bin5line->Draw("SAME");
    TLegend *legend_mu_loose_1d=new TLegend(.7,.7,.85,.85);
    legend_mu_loose_1d->AddEntry(syst_mu_sf_loose_eta1,"| #eta | < 1.479", "lp");
    legend_mu_loose_1d->AddEntry(syst_mu_sf_loose_eta2,"| #eta | > 1.479", "lp");
    legend_mu_loose_1d->SetFillColor(0);
    legend_mu_loose_1d->Draw("SAME");
    c_syst_mu_sf_loose_1d->Update();
    c_syst_mu_sf_loose_1d->Print((plots_dir+"syst_mu_sf_loose_1d.png").c_str());
 
    TCanvas *c_stat_mu_sf_loose = new TCanvas("c_stat_mu_sf_loose", "Stat. unc. on mu. loose SF", 400, 600);
    stat_mu_sf_loose->Draw("TEXT COLZ");
    gPad->Update();
    stat_mu_sf_loose->SetTitle("Stat. unc. on mu. loose SF");
    stat_mu_sf_loose->GetXaxis()->SetTitle("| #eta |");
    stat_mu_sf_loose->GetXaxis()->SetTitleOffset(0.9);
    stat_mu_sf_loose->GetXaxis()->SetTitleSize(0.04);
    stat_mu_sf_loose->GetXaxis()->SetLabelSize(0.02);
    stat_mu_sf_loose->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_mu_sf_loose->GetYaxis()->SetTitleOffset(0.9);
    stat_mu_sf_loose->GetYaxis()->SetTitleSize(0.04);
    stat_mu_sf_loose->GetYaxis()->SetLabelSize(0.02);
    stat_mu_sf_loose->GetYaxis()->SetRangeUser(10,100);
    stat_mu_sf_loose->SetMinimum(0);
    stat_mu_sf_loose->SetMaximum(0.2);
    stat_mu_sf_loose->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_mu_sf_loose->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_mu_sf_loose->Update();
    c_stat_mu_sf_loose->Print((plots_dir+"stat_mu_sf_loose.png").c_str());

    // mu medium plots
    TCanvas *c_syst_mu_sf_signal_medium = new TCanvas("c_syst_mu_sf_signal_medium", "Syst. unc. on mu. medium SF (signal)", 400, 600);
    syst_mu_sf_signal_medium->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_signal_medium->SetTitle("Syst. unc. on mu. medium SF from signal choice");
    syst_mu_sf_signal_medium->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_signal_medium->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_medium->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_medium->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_medium->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_signal_medium->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_medium->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_medium->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_medium->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_signal_medium->SetMinimum(0);
    syst_mu_sf_signal_medium->SetMaximum(0.2);
    syst_mu_sf_signal_medium->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_signal_medium->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_signal_medium->Update();
    c_syst_mu_sf_signal_medium->Print((plots_dir+"syst_mu_sf_signal_medium.png").c_str());
    
    TCanvas *c_syst_mu_sf_background_medium = new TCanvas("c_syst_mu_sf_background_medium", "Syst. unc. on mu. medium SF (b.g.)", 400, 600);
    syst_mu_sf_background_medium->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_background_medium->SetTitle("Syst. unc. on mu. medium SF from background choice");
    syst_mu_sf_background_medium->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_background_medium->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_medium->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_medium->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_medium->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_background_medium->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_medium->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_medium->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_medium->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_background_medium->SetMinimum(0);
    syst_mu_sf_background_medium->SetMaximum(0.2);
    syst_mu_sf_background_medium->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_background_medium->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_background_medium->Update();
    c_syst_mu_sf_background_medium->Print((plots_dir+"syst_mu_sf_background_medium.png").c_str());
    
    TCanvas *c_syst_mu_sf_combined_medium = new TCanvas("c_syst_mu_sf_combined_medium", "Syst. unc. on mu. medium SF (combined)", 400, 600);
    syst_mu_sf_combined_medium->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_combined_medium->SetTitle("Syst. unc. on mu. medium SF (combined)");
    syst_mu_sf_combined_medium->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_combined_medium->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_medium->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_medium->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_medium->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_combined_medium->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_medium->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_medium->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_medium->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_combined_medium->SetMinimum(0);
    syst_mu_sf_combined_medium->SetMaximum(0.2);
    syst_mu_sf_combined_medium->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_combined_medium->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_combined_medium->Update();
    c_syst_mu_sf_combined_medium->Print((plots_dir+"syst_mu_sf_combined_medium.png").c_str());
    
    TCanvas *c_syst_mu_sf_medium_1d = new TCanvas("c_syst_mu_sf_medium_1d","Systematic uncertainty on medium muon scale factors",800,600);
    c_syst_mu_sf_medium_1d->SetGrid(0,1);
    TH1D *syst_mu_sf_medium_eta1 = syst_mu_sf_combined_medium->ProjectionY("syst_mu_sf_medium_eta1", 1, 1);
    TH1D *syst_mu_sf_medium_eta2 = syst_mu_sf_combined_medium->ProjectionY("syst_mu_sf_medium_eta2", 2, 2);
    syst_mu_sf_medium_eta1->Draw("P0 L");
    syst_mu_sf_medium_eta1->SetTitle("Systematic uncertainty on medium muon scale factors");
    syst_mu_sf_medium_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_medium_eta1->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_medium_eta1->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_medium_eta1->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_medium_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
    syst_mu_sf_medium_eta1->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_medium_eta1->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_medium_eta1->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_medium_eta1->SetMinimum(0);
    syst_mu_sf_medium_eta1->SetMaximum(0.1);
    syst_mu_sf_medium_eta1->SetMarkerStyle(20);
    syst_mu_sf_medium_eta1->SetMarkerColor(1);
    syst_mu_sf_medium_eta1->SetLineColor(1);
    syst_mu_sf_medium_eta2->SetMarkerStyle(21);
    syst_mu_sf_medium_eta2->SetMarkerColor(1861);
    syst_mu_sf_medium_eta2->SetLineColor(1861);
    syst_mu_sf_medium_eta2->Draw("P0 L SAME");
    bin1line->Draw("SAME");
    bin2line->Draw("SAME");
    bin3line->Draw("SAME");
    bin4line->Draw("SAME");
    bin5line->Draw("SAME");
    TLegend *legend_mu_medium_1d=new TLegend(.7,.7,.85,.85);
    legend_mu_medium_1d->AddEntry(syst_mu_sf_medium_eta1,"| #eta | < 1.479", "lp");
    legend_mu_medium_1d->AddEntry(syst_mu_sf_medium_eta2,"| #eta | > 1.479", "lp");
    legend_mu_medium_1d->SetFillColor(0);
    legend_mu_medium_1d->Draw("SAME");
    c_syst_mu_sf_medium_1d->Update();
    c_syst_mu_sf_medium_1d->Print((plots_dir+"syst_mu_sf_medium_1d.png").c_str());

    TCanvas *c_stat_mu_sf_medium = new TCanvas("c_stat_mu_sf_medium", "Stat. unc. on mu. medium SF", 400, 600);
    stat_mu_sf_medium->Draw("TEXT COLZ");
    gPad->Update();
    stat_mu_sf_medium->SetTitle("Stat. unc. on mu. medium SF");
    stat_mu_sf_medium->GetXaxis()->SetTitle("| #eta |");
    stat_mu_sf_medium->GetXaxis()->SetTitleOffset(0.9);
    stat_mu_sf_medium->GetXaxis()->SetTitleSize(0.04);
    stat_mu_sf_medium->GetXaxis()->SetLabelSize(0.02);
    stat_mu_sf_medium->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_mu_sf_medium->GetYaxis()->SetTitleOffset(0.9);
    stat_mu_sf_medium->GetYaxis()->SetTitleSize(0.04);
    stat_mu_sf_medium->GetYaxis()->SetLabelSize(0.02);
    stat_mu_sf_medium->GetYaxis()->SetRangeUser(10,100);
    stat_mu_sf_medium->SetMinimum(0);
    stat_mu_sf_medium->SetMaximum(0.2);
    stat_mu_sf_medium->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_mu_sf_medium->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_mu_sf_medium->Update();
    c_stat_mu_sf_medium->Print((plots_dir+"stat_mu_sf_medium.png").c_str());

    // mu tight plots
    TCanvas *c_syst_mu_sf_signal_tight = new TCanvas("c_syst_mu_sf_signal_tight", "Syst. unc. on mu. tight SF (signal)", 400, 600);
    syst_mu_sf_signal_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_signal_tight->SetTitle("Syst. unc. on mu. tight SF from signal choice");
    syst_mu_sf_signal_tight->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_signal_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_tight->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_tight->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_signal_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_signal_tight->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_signal_tight->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_signal_tight->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_signal_tight->SetMinimum(0);
    syst_mu_sf_signal_tight->SetMaximum(0.2);
    syst_mu_sf_signal_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_signal_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_signal_tight->Update();
    c_syst_mu_sf_signal_tight->Print((plots_dir+"syst_mu_sf_signal_tight.png").c_str());
    
    TCanvas *c_syst_mu_sf_background_tight = new TCanvas("c_syst_mu_sf_background_tight", "Syst. unc. on mu. tight SF (b.g.)", 400, 600);
    syst_mu_sf_background_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_background_tight->SetTitle("Syst. unc. on mu. tight SF from background choice");
    syst_mu_sf_background_tight->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_background_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_tight->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_tight->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_background_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_background_tight->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_background_tight->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_background_tight->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_background_tight->SetMinimum(0);
    syst_mu_sf_background_tight->SetMaximum(0.2);
    syst_mu_sf_background_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_background_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_background_tight->Update();
    c_syst_mu_sf_background_tight->Print((plots_dir+"syst_mu_sf_background_tight.png").c_str());
    
    TCanvas *c_syst_mu_sf_combined_tight = new TCanvas("c_syst_mu_sf_combined_tight", "Syst. unc. on mu. tight SF (combined)", 400, 600);
    syst_mu_sf_combined_tight->Draw("TEXT COLZ");
    gPad->Update();
    syst_mu_sf_combined_tight->SetTitle("Syst. unc. on mu. tight SF (combined)");
    syst_mu_sf_combined_tight->GetXaxis()->SetTitle("| #eta |");
    syst_mu_sf_combined_tight->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_tight->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_tight->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_combined_tight->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_combined_tight->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_combined_tight->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_combined_tight->GetYaxis()->SetRangeUser(10,100);
    syst_mu_sf_combined_tight->SetMinimum(0);
    syst_mu_sf_combined_tight->SetMaximum(0.2);
    syst_mu_sf_combined_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_mu_sf_combined_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_mu_sf_combined_tight->Update();
    c_syst_mu_sf_combined_tight->Print((plots_dir+"syst_mu_sf_combined_tight.png").c_str());
    
    TCanvas *c_syst_mu_sf_tight_1d = new TCanvas("c_syst_mu_sf_tight_1d","Systematic uncertainty on tight muon scale factors",800,600);
    c_syst_mu_sf_tight_1d->SetGrid(0,1);
    TH1D *syst_mu_sf_tight_eta1 = syst_mu_sf_combined_tight->ProjectionY("syst_mu_sf_tight_eta1", 1, 1);
    TH1D *syst_mu_sf_tight_eta2 = syst_mu_sf_combined_tight->ProjectionY("syst_mu_sf_tight_eta2", 2, 2);
    syst_mu_sf_tight_eta1->Draw("P0 L");
    syst_mu_sf_tight_eta1->SetTitle("Systematic uncertainty on tight muon scale factors");
    syst_mu_sf_tight_eta1->GetXaxis()->SetTitle("p_{T} [GeV]");
    syst_mu_sf_tight_eta1->GetXaxis()->SetTitleOffset(0.9);
    syst_mu_sf_tight_eta1->GetXaxis()->SetTitleSize(0.04);
    syst_mu_sf_tight_eta1->GetXaxis()->SetLabelSize(0.02);
    syst_mu_sf_tight_eta1->GetYaxis()->SetTitle("#varepsilon_{data} / #varepsilon_{MC}");
    syst_mu_sf_tight_eta1->GetYaxis()->SetTitleOffset(0.9);
    syst_mu_sf_tight_eta1->GetYaxis()->SetTitleSize(0.04);
    syst_mu_sf_tight_eta1->GetYaxis()->SetLabelSize(0.02);
    syst_mu_sf_tight_eta1->SetMinimum(0);
    syst_mu_sf_tight_eta1->SetMaximum(0.1);
    syst_mu_sf_tight_eta1->SetMarkerStyle(20);
    syst_mu_sf_tight_eta1->SetMarkerColor(1);
    syst_mu_sf_tight_eta1->SetLineColor(1);
    syst_mu_sf_tight_eta2->SetMarkerStyle(21);
    syst_mu_sf_tight_eta2->SetMarkerColor(1861);
    syst_mu_sf_tight_eta2->SetLineColor(1861);
    syst_mu_sf_tight_eta2->Draw("P0 L SAME");
    bin1line->Draw("SAME");
    bin2line->Draw("SAME");
    bin3line->Draw("SAME");
    bin4line->Draw("SAME");
    bin5line->Draw("SAME");
    TLegend *legend_mu_tight_1d=new TLegend(.7,.7,.85,.85);
    legend_mu_tight_1d->AddEntry(syst_mu_sf_tight_eta1,"| #eta | < 1.479", "lp");
    legend_mu_tight_1d->AddEntry(syst_mu_sf_tight_eta2,"| #eta | > 1.479", "lp");
    legend_mu_tight_1d->SetFillColor(0);
    legend_mu_tight_1d->Draw("SAME");
    c_syst_mu_sf_tight_1d->Update();
    c_syst_mu_sf_tight_1d->Print((plots_dir+"syst_mu_sf_tight_1d.png").c_str());

    TCanvas *c_stat_mu_sf_tight = new TCanvas("c_stat_mu_sf_tight", "Stat. unc. on mu. tight SF", 400, 600);
    stat_mu_sf_tight->Draw("TEXT COLZ");
    gPad->Update();
    stat_mu_sf_tight->SetTitle("Stat. unc. on mu. tight SF");
    stat_mu_sf_tight->GetXaxis()->SetTitle("| #eta |");
    stat_mu_sf_tight->GetXaxis()->SetTitleOffset(0.9);
    stat_mu_sf_tight->GetXaxis()->SetTitleSize(0.04);
    stat_mu_sf_tight->GetXaxis()->SetLabelSize(0.02);
    stat_mu_sf_tight->GetYaxis()->SetTitle("p_{T} [GeV]");
    stat_mu_sf_tight->GetYaxis()->SetTitleOffset(0.9);
    stat_mu_sf_tight->GetYaxis()->SetTitleSize(0.04);
    stat_mu_sf_tight->GetYaxis()->SetLabelSize(0.02);
    stat_mu_sf_tight->GetYaxis()->SetRangeUser(10,100);
    stat_mu_sf_tight->SetMinimum(0);
    stat_mu_sf_tight->SetMaximum(0.2);
    stat_mu_sf_tight->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) stat_mu_sf_tight->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_stat_mu_sf_tight->Update();
    c_stat_mu_sf_tight->Print((plots_dir+"stat_mu_sf_tight.png").c_str());


  }
  
  TFile *uncertainties_ele_veto = new TFile((root_dir+"combined_uncertainties_ele_veto.root").c_str(),"RECREATE");
  syst_ele_sf_signal_veto     ->Write("syst_ele_sf_signal_veto");
  syst_ele_sf_background_veto ->Write("syst_ele_sf_background_veto");
  syst_ele_sf_combined_veto   ->Write("syst_ele_sf_combined_veto");
  stat_ele_sf_veto            ->Write("stat_ele_sf_veto");
  TFile *uncertainties_ele_loose = new TFile((root_dir+"combined_uncertainties_ele_loose.root").c_str(),"RECREATE");
  syst_ele_sf_signal_loose     ->Write("syst_ele_sf_signal_loose");
  syst_ele_sf_background_loose ->Write("syst_ele_sf_background_loose");
  syst_ele_sf_combined_loose   ->Write("syst_ele_sf_combined_loose");
  stat_ele_sf_loose            ->Write("stat_ele_sf_loose");
  TFile *uncertainties_ele_medium = new TFile((root_dir+"combined_uncertainties_ele_medium.root").c_str(),"RECREATE");
  syst_ele_sf_signal_medium     ->Write("syst_ele_sf_signal_medium");
  syst_ele_sf_background_medium ->Write("syst_ele_sf_background_medium");
  syst_ele_sf_combined_medium   ->Write("syst_ele_sf_combined_medium");
  stat_ele_sf_medium            ->Write("stat_ele_sf_medium");
  TFile *uncertainties_ele_tight = new TFile((root_dir+"combined_uncertainties_ele_tight.root").c_str(),"RECREATE");
  syst_ele_sf_signal_tight     ->Write("syst_ele_sf_signal_tight");
  syst_ele_sf_background_tight ->Write("syst_ele_sf_background_tight");
  syst_ele_sf_combined_tight   ->Write("syst_ele_sf_combined_tight");
  stat_ele_sf_tight            ->Write("stat_ele_sf_tight");
  TFile *uncertainties_mu_veto = new TFile((root_dir+"combined_uncertainties_mu_veto.root").c_str(),"RECREATE");
  syst_mu_sf_signal_veto     ->Write("syst_mu_sf_signal_veto");
  syst_mu_sf_background_veto ->Write("syst_mu_sf_background_veto");
  syst_mu_sf_combined_veto   ->Write("syst_mu_sf_combined_veto");
  stat_mu_sf_veto            ->Write("stat_mu_sf_veto");
  TFile *uncertainties_mu_loose = new TFile((root_dir+"combined_uncertainties_mu_loose.root").c_str(),"RECREATE");
  syst_mu_sf_signal_loose     ->Write("syst_mu_sf_signal_loose");
  syst_mu_sf_background_loose ->Write("syst_mu_sf_background_loose");
  syst_mu_sf_combined_loose   ->Write("syst_mu_sf_combined_loose");
  stat_mu_sf_loose            ->Write("stat_mu_sf_loose");
  TFile *uncertainties_mu_medium = new TFile((root_dir+"combined_uncertainties_mu_medium.root").c_str(),"RECREATE");
  syst_mu_sf_signal_medium     ->Write("syst_mu_sf_signal_medium");
  syst_mu_sf_background_medium ->Write("syst_mu_sf_background_medium");
  syst_mu_sf_combined_medium   ->Write("syst_mu_sf_combined_medium");
  stat_mu_sf_medium            ->Write("stat_mu_sf_medium");
  TFile *uncertainties_mu_tight = new TFile((root_dir+"combined_uncertainties_mu_tight.root").c_str(),"RECREATE");
  syst_mu_sf_signal_tight     ->Write("syst_mu_sf_signal_tight");
  syst_mu_sf_background_tight ->Write("syst_mu_sf_background_tight");
  syst_mu_sf_combined_tight   ->Write("syst_mu_sf_combined_tight");
  stat_mu_sf_tight            ->Write("stat_mu_sf_tight");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void flavor_yield_ratio(
  string plots_dir,
  string root_dir,
  bool verbose=false
) {
  gStyle->SetOptStat(0); 
  int n_skims=4; 
  int nbins=17;
  Float_t Z_pT_bins[] = {0., 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 300, 400};
  int max_pT = 400;
  TFile *f_sf_ele  = TFile::Open((root_dir+"scalefactors_ele_74x.root").c_str(),"READ");
  TFile *f_sf_mu   = TFile::Open((root_dir+"scalefactors_mu_74x.root").c_str(),"READ");
  //TFile *f_uncertainties_ele_tight  = TFile::Open((root_dir+"combined_uncertainties_ele_tight.root").c_str(),"READ");
  //TFile *f_uncertainties_mu_tight   = TFile::Open((root_dir+"combined_uncertainties_mu_tight.root").c_str(),"READ");
  TFile *f_mu_triggers = TFile::Open("/home/dhsu/TagAndProbe/Data_Triggers/IsoMu27/eff.root", "READ");
  TFile *f_ele_triggers = TFile::Open("/home/dhsu/TagAndProbe/Data_Triggers/Ele27_eta2p1_WPLoose_Gsf/eff.root", "READ");
  
  TClonesArray skim_files ("TFile", n_skims);
  new(skim_files[0]) TFile("~/leptonScaleFactors/root/DYJetsToLL_74x_aMCatNLO_genMatching_BaselineToTight_electronTnP.root","READ");
  new(skim_files[1]) TFile("~/leptonScaleFactors/root/DYJetsToLL_74x_aMCatNLO_genMatching_BaselineToTight_muonTnP.root","READ");
  new(skim_files[2]) TFile("~/leptonScaleFactors/root/SingleElectron_BaselineToTight_electronTnP.root","READ");
  new(skim_files[3]) TFile("~/leptonScaleFactors/root/SingleMuon_BaselineToTight_muonTnP.root","READ");

  TH2D *syst_ele_sf_tight    = (TH2D*) f_sf_ele->Get("scalefactors_Tight_ele_syst_error_combined");
  TH2D *stat_ele_sf_tight_lo = (TH2D*) f_sf_ele->Get("scalefactors_Tight_ele_stat_error_lo");
  TH2D *stat_ele_sf_tight_hi = (TH2D*) f_sf_ele->Get("scalefactors_Tight_ele_stat_error_hi");
  TH2D *syst_mu_sf_tight    = (TH2D*) f_sf_mu->Get("scalefactors_Tight_mu_syst_error_combined");
  TH2D *stat_mu_sf_tight_lo = (TH2D*) f_sf_mu->Get("scalefactors_Tight_mu_stat_error_lo");
  TH2D *stat_mu_sf_tight_hi = (TH2D*) f_sf_mu->Get("scalefactors_Tight_mu_stat_error_hi");

  TH2D *trig_eff_ele = (TH2D*) f_ele_triggers->Get("hEffEtaPt");
  TH2D *trig_eff_mu  = (TH2D*) f_mu_triggers->Get("hEffEtaPt");
  TH2D *unc_trig_eff_ele = (TH2D*) f_ele_triggers->Get("hErrhEtaPt");
  TH2D *unc_trig_eff_mu  = (TH2D*) f_mu_triggers->Get("hErrhEtaPt");
  TH2D *sf_ele_tight = (TH2D*) f_sf_ele->Get("scalefactors_Tight_ele");
  TH2D *sf_mu_tight  = (TH2D*) f_sf_mu->Get("scalefactors_Tight_mu");
  //printf("%f %f\n", sf_ele_tight->GetBinContent(sf_ele_tight->FindBin(1.3,45)), sf_mu_tight->GetRMS());
  //return;
  TClonesArray dN_dZpT_ ("TH1D", n_skims);
  new(dN_dZpT_[0]) TH1D("ZpT_mc_dN_dZpT_ele_tight",   "p_{T}^{ll} distribution of Z #rightarrow ee events in MC", nbins, Z_pT_bins);
  new(dN_dZpT_[1]) TH1D("ZpT_mc_dN_dZpT_mu_tight",    "p_{T}^{ll} distribution of Z #rightarrow #mu#mu events in MC", nbins, Z_pT_bins);
  new(dN_dZpT_[2]) TH1D("ZpT_data_dN_dZpT_ele_tight", "p_{T}^{ll} distribution of Z #rightarrow ee events in data", nbins, Z_pT_bins);
  new(dN_dZpT_[3]) TH1D("ZpT_data_dN_dZpT_mu_tight",  "p_{T}^{ll} distribution of Z #rightarrow #mu#mu events in data", nbins, Z_pT_bins);
  TClonesArray syst_sf_ZpT_ ("TH1D", n_skims);
  new(syst_sf_ZpT_[0]) TH1D("ZpT_mc_syst_sf_ele_tight",   "S.F. systematics as p_{T} of Z #rightarrow ee events in MC", nbins, Z_pT_bins);
  new(syst_sf_ZpT_[1]) TH1D("ZpT_mc_syst_sf_mu_tight",    "S.F. systematics as p_{T} of Z #rightarrow #mu#mu events in MC", nbins, Z_pT_bins);
  new(syst_sf_ZpT_[2]) TH1D("ZpT_data_syst_sf_ele_tight", "S.F. systematics as p_{T} of Z #rightarrow ee events in data", nbins, Z_pT_bins);
  new(syst_sf_ZpT_[3]) TH1D("ZpT_data_syst_sf_mu_tight",  "S.F. systematics as p_{T} of Z #rightarrow #mu#mu events in data", nbins, Z_pT_bins);
  TClonesArray stat_sf_ZpT_ ("TH1D", n_skims);
  new(stat_sf_ZpT_[0]) TH1D("ZpT_mc_stat_sf_ele_tight",   "S.F. stat. unc. as p_{T} of Z #rightarrow ee events in MC", nbins, Z_pT_bins);
  new(stat_sf_ZpT_[1]) TH1D("ZpT_mc_stat_sf_mu_tight",    "S.F. stat. unc. as p_{T} of Z #rightarrow #mu#mu events in MC", nbins, Z_pT_bins);
  new(stat_sf_ZpT_[2]) TH1D("ZpT_data_stat_sf_ele_tight", "S.F. stat. unc. as p_{T} of Z #rightarrow ee events in data", nbins, Z_pT_bins);
  new(stat_sf_ZpT_[3]) TH1D("ZpT_data_stat_sf_mu_tight",  "S.F. stat. unc. as p_{T} of Z #rightarrow #mu#mu events in data", nbins, Z_pT_bins);
  TClonesArray trig_unc_ZpT_ ("TH1D", n_skims);
  new(trig_unc_ZpT_[0]) TH1D("ZpT_mc_trig_unc_ele_tight",   "Trigger uncertainty as p_{T} of Z #rightarrow ee events in MC", nbins, Z_pT_bins);
  new(trig_unc_ZpT_[1]) TH1D("ZpT_mc_trig_unc_mu_tight",    "Trigger uncertainty as p_{T} of Z #rightarrow #mu#mu events in MC", nbins, Z_pT_bins);
  new(trig_unc_ZpT_[2]) TH1D("ZpT_data_trig_unc_ele_tight", "Trigger uncertainty as p_{T} of Z #rightarrow ee events in data", nbins, Z_pT_bins);
  new(trig_unc_ZpT_[3]) TH1D("ZpT_data_trig_unc_mu_tight",  "Trigger uncertainty as p_{T} of Z #rightarrow #mu#mu events in data", nbins, Z_pT_bins);
  TClonesArray stat_events_ZpT_ ("TH1D", n_skims);
  new(stat_events_ZpT_[0]) TH1D("ZpT_mc_stat_events_ele_tight",   "Stat. unc. as p_{T} of Z #rightarrow ee events in MC", nbins, Z_pT_bins);
  new(stat_events_ZpT_[1]) TH1D("ZpT_mc_stat_events_mu_tight",    "Stat. unc. as p_{T} of Z #rightarrow #mu#mu events in MC", nbins, Z_pT_bins);
  new(stat_events_ZpT_[2]) TH1D("ZpT_data_stat_events_ele_tight", "Stat. unc. as p_{T} of Z #rightarrow ee events in data", nbins, Z_pT_bins);
  new(stat_events_ZpT_[3]) TH1D("ZpT_data_stat_events_mu_tight",  "Stat. unc. as p_{T} of Z #rightarrow #mu#mu events in data", nbins, Z_pT_bins);
  TClonesArray unc_ZpT_ ("TH1D", n_skims);
  new(unc_ZpT_[0]) TH1D("ZpT_mc_unc_ele_tight",   "Uncertainty as p_{T} of Z #rightarrow ee events in MC", nbins, Z_pT_bins);
  new(unc_ZpT_[1]) TH1D("ZpT_mc_unc_mu_tight",    "Uncertainty as p_{T} of Z #rightarrow #mu#mu events in MC", nbins, Z_pT_bins);
  new(unc_ZpT_[2]) TH1D("ZpT_data_unc_ele_tight", "Uncertainty as p_{T} of Z #rightarrow ee events in data", nbins, Z_pT_bins);
  new(unc_ZpT_[3]) TH1D("ZpT_data_unc_mu_tight",  "Uncertainty as p_{T} of Z #rightarrow #mu#mu events in data", nbins, Z_pT_bins);
  //TClonesArray sum_weights_vs_ZpT_ ("TH1D", n_skims);
  //new(sum_weights_vs_ZpT_[0]) TH1D("h_mc_sum_weights_ele_tight",   "p_{T} distribution of Z #rightarrow ee events in MC", nbins, Z_pT_bins);
  //new(sum_weights_vs_ZpT_[1]) TH1D("h_mc_sum_weights_mu_tight",    "p_{T} distribution of Z #rightarrow {#mu}{#mu} events in MC", nbins, Z_pT_bins);
  //new(sum_weights_vs_ZpT_[2]) TH1D("h_data_sum_weights_ele_tight", "p_{T} distribution of Z #rightarrow ee events in data", nbins, Z_pT_bins);
  //new(sum_weights_vs_ZpT_[3]) TH1D("h_data_sum_weights_mu_tight",  "p_{T} distribution of Z #rightarrow {#mu}{#mu} events in data", nbins, Z_pT_bins);
  
  // read from tnp skim
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
  
  for(int skim=0; skim<n_skims; skim++) {
    printf("looping over skim %d . . .\n", skim);
    TTree *tree = (TTree*) ((TFile*)skim_files[skim])->Get("Events");
    tree->SetBranchAddress("runNum",   &runNum   );  
    tree->SetBranchAddress("lumiSec",  &lumiSec  );  
    tree->SetBranchAddress("evtNum",   &evtNum   );  
    tree->SetBranchAddress("npv",      &npv      );  
    tree->SetBranchAddress("pass",     &pass     );  
    tree->SetBranchAddress("npu",      &npu      );  
    tree->SetBranchAddress("scale1fb", &scale1fb );
    tree->SetBranchAddress("mass",     &mass     );  
    tree->SetBranchAddress("qtag",     &qtag     );  
    tree->SetBranchAddress("qprobe",   &qprobe   );  
    tree->SetBranchAddress("tag",      &p4_tag   );  
    tree->SetBranchAddress("probe",    &p4_probe );     
    Long64_t nentries= tree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {
      tree->GetEntry(i);
      if(!(
        pass==1 &&
        p4_tag->Pt() >= 30 &&
        TMath::Abs(p4_tag->Eta()) <= 2.1 &&
        TMath::Abs(mass - 90) <= 30 &&
        qtag + qprobe == 0
      )) continue;
      TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
      double Z_pT = systemP4.Pt();
      double trig_eff_tag=1, trig_eff_probe=1, sf_tag=1, sf_probe=1;
      double stat_sf_tag=0, stat_sf_probe=0, syst_sf_tag=0, syst_sf_probe=0, unc_trig_tag=0, unc_trig_probe;
      double tag_eta = p4_tag->Eta();
      double tag_pT = p4_tag->Pt();
      double probe_eta = p4_probe->Eta();
      double probe_pT = p4_probe->Pt();
      double mc_reweighting=1;
      if(tag_eta >= 2.5)    continue;
      if(tag_eta <= -2.5)   continue;
      if(probe_eta >= 2.5)  continue;
      if(probe_eta <= -2.5) continue;
      //if(tag_pT >= 200) tag_pT = 199.9;
      //if(probe_pT >= 200) probe_pT = 199.9;
      //if(tag_eta >= 2.5) tag_eta = 2.49;
      //if(tag_eta <= -2.5) tag_eta = -2.49;
      //if(probe_eta >= 2.5) probe_eta = 2.49;
      //if(probe_eta <= -2.5) probe_eta = -2.49;
      switch(skim) {
        case 0:
          if(probe_pT >= 200) probe_pT = 199.9;
          if(tag_pT >= 200) tag_pT = 199.9;
          sf_tag = sf_ele_tight->GetBinContent( sf_ele_tight->FindBin( TMath::Abs(tag_eta), tag_pT ));
          sf_probe = sf_ele_tight->GetBinContent( sf_ele_tight->FindBin( TMath::Abs(probe_eta), probe_pT ));
          
          stat_sf_tag   = TMath::Max(
            stat_ele_sf_tight_hi->GetBinContent(stat_ele_sf_tight_hi->FindBin( TMath::Abs(tag_eta), tag_pT)), 
            stat_ele_sf_tight_lo->GetBinContent(stat_ele_sf_tight_lo->FindBin( TMath::Abs(tag_eta), tag_pT)) 
          );
          stat_sf_probe = TMath::Max(
            stat_ele_sf_tight_hi->GetBinContent(stat_ele_sf_tight_hi->FindBin( TMath::Abs(probe_eta), probe_pT)),
            stat_ele_sf_tight_lo->GetBinContent(stat_ele_sf_tight_lo->FindBin( TMath::Abs(probe_eta), probe_pT))
          );
          syst_sf_tag   = syst_ele_sf_tight->GetBinContent(syst_ele_sf_tight->FindBin( TMath::Abs(tag_eta), tag_pT));
          syst_sf_probe = syst_ele_sf_tight->GetBinContent(syst_ele_sf_tight->FindBin( TMath::Abs(probe_eta), probe_pT));
          
          trig_eff_tag   = trig_eff_ele->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(tag_eta), tag_pT ) );
          if(probe_pT < 30 || TMath::Abs(probe_eta) > 2.1) trig_eff_probe = 0;
          else              trig_eff_probe = trig_eff_ele->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(probe_eta), probe_pT ) );

          unc_trig_tag   = unc_trig_eff_ele->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(tag_eta),  tag_pT ) );
          if(probe_pT < 30 || TMath::Abs(probe_eta) > 2.1) unc_trig_probe = 0;
          else              unc_trig_probe = unc_trig_eff_ele->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(probe_eta), probe_pT ) );
          mc_reweighting = 2151*6025.2/19312436.;
          break;
        case 1:
          if(tag_eta >= 2.4) tag_eta = 2.39;
          if(tag_eta <= -2.4) tag_eta = -2.39;
          if(probe_eta >= 2.4) probe_eta = 2.39;
          if(probe_eta <= -2.4) probe_eta = -2.39;
          sf_tag = sf_mu_tight->GetBinContent( sf_mu_tight->FindBin( TMath::Abs(tag_eta), tag_pT ));
          sf_probe = sf_mu_tight->GetBinContent( sf_mu_tight->FindBin( TMath::Abs(probe_eta), probe_pT ));
          stat_sf_tag   = TMath::Max(
            stat_mu_sf_tight_hi->GetBinContent(stat_mu_sf_tight_hi->FindBin( TMath::Abs(tag_eta), tag_pT)), 
            stat_mu_sf_tight_lo->GetBinContent(stat_mu_sf_tight_lo->FindBin( TMath::Abs(tag_eta), tag_pT)) 
          );
          stat_sf_probe = TMath::Max(
            stat_mu_sf_tight_hi->GetBinContent(stat_mu_sf_tight_hi->FindBin( TMath::Abs(probe_eta), probe_pT)),
            stat_mu_sf_tight_lo->GetBinContent(stat_mu_sf_tight_lo->FindBin( TMath::Abs(probe_eta), probe_pT))
          );
          syst_sf_tag   = syst_mu_sf_tight->GetBinContent(syst_mu_sf_tight->FindBin( TMath::Abs(tag_eta), tag_pT));
          syst_sf_probe = syst_mu_sf_tight->GetBinContent(syst_mu_sf_tight->FindBin( TMath::Abs(probe_eta), probe_pT));

          trig_eff_tag   = trig_eff_mu->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(tag_eta), tag_pT) );
          if(probe_pT < 30 || TMath::Abs(probe_eta) > 2.1) trig_eff_probe = 0;
          else              trig_eff_probe = trig_eff_mu->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(probe_eta), probe_pT) );

          unc_trig_tag   = unc_trig_eff_mu->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(tag_eta), tag_pT) );
          if(probe_pT < 30 || TMath::Abs(probe_eta) > 2.1) unc_trig_probe = 0;
          else              unc_trig_probe = unc_trig_eff_mu->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(probe_eta), probe_pT) );
          mc_reweighting = 2106*6025.2/19312436.;
          break;
        case 2:
          break;
        case 3:
          break;
        default:
          break;
      }
      //if(trig_eff_tag==0 || !(unc_trig_tag < 1) || !(unc_trig_tag>0)) continue;
      // scale1fb encapsulates pileup reweighting and generator weights
      double weight = sf_tag * sf_probe * (1 - (1 - trig_eff_tag)*(1 - trig_eff_probe)) * scale1fb * mc_reweighting;
      if(skim>1) weight=1; //lumi
      if(trig_eff_tag==0) printf("tag trig eff = 0 for tag (%f, %f)  probe (%f, %f)\n", p4_tag->Eta(), p4_tag->Pt(), p4_probe->Eta(), p4_probe->Pt());
      //if(trig_eff_probe==0) printf("probe trig eff = 0 for tag (%f, %f)  probe (%f, %f)\n", p4_tag->Eta(), p4_tag->Pt(), p4_probe->Eta(), p4_probe->Pt());
      if(sf_tag == 0) printf("sf_tag = 0 for tag (%f, %f)  probe (%f, %f)\n", p4_tag->Eta(), p4_tag->Pt(), p4_probe->Eta(), p4_probe->Pt());
      if(sf_probe == 0) printf("sf_probe = 0 for tag (%f, %f)  probe (%f, %f)\n", p4_tag->Eta(), p4_tag->Pt(), p4_probe->Eta(), p4_probe->Pt());
      // syst. unc. for tag and probe are 100% correlated
      // stat. unc. are not
      ( (TH1D*) dN_dZpT_[skim]    )->Fill(Z_pT, weight);
      if(verbose && skim<2) printf("filling Z_pT %f with N %f , sum is now %f \n", Z_pT, weight, ((TH1D*)dN_dZpT_[skim])->GetBinContent(((TH1D*)dN_dZpT_[skim])->FindBin(Z_pT)));
      ( (TH1D*) syst_sf_ZpT_[skim] )->Fill(Z_pT, weight*( syst_sf_tag + syst_sf_probe) / (sf_tag * sf_probe) );
      if(verbose && skim<2) printf("filling Z_pT %f with weighted syst %f , sum is now %f \n", Z_pT, weight*(syst_sf_tag+syst_sf_probe)/(sf_tag*sf_probe), ((TH1D*)syst_sf_ZpT_[skim])->GetBinContent(((TH1D*)syst_sf_ZpT_[skim])->FindBin(Z_pT)));
      // take square root of quadrature sum of stat. unc. after loop
      ( (TH1D*) stat_events_ZpT_[skim] )->Fill(Z_pT, weight*weight);
      ( (TH1D*) stat_sf_ZpT_[skim] )->Fill(Z_pT, weight*weight*( pow(stat_sf_tag/sf_tag,2) + pow(stat_sf_probe/sf_probe,2) ));
      
      // trigger uncertainties
      if(trig_eff_probe!=0) ( (TH1D*) trig_unc_ZpT_[skim] )->Fill(Z_pT, weight*weight*( pow(unc_trig_tag/trig_eff_tag, 2) + pow(unc_trig_probe/trig_eff_probe, 2)));
      else                  ( (TH1D*) trig_unc_ZpT_[skim] )->Fill(Z_pT, pow(weight* unc_trig_tag / trig_eff_tag,2));
    }
  }
  //double ele_scale = ((TH1D*)dN_dZpT_[2])->Integral() /  ((TH1D*)dN_dZpT_[0])->Integral();
  //double mu_scale  = ((TH1D*)dN_dZpT_[3])->Integral() /  ((TH1D*)dN_dZpT_[1])->Integral();
  //printf("Scaling ele MC by factor of %f\n", ele_scale);
  //printf("Scaling mu MC by factor of %f\n", mu_scale);

  //((TH1D*)dN_dZpT_[0])->Scale(ele_scale);
  //((TH1D*)dN_dZpT_[1])->Scale(mu_scale);
  
  for(int skim=0; skim<n_skims; skim++) {
    for(int j=1; j<=nbins; j++) {
      ((TH1D*)unc_ZpT_[skim])->SetBinContent(j,
        sqrt(
          pow(((TH1D*)syst_sf_ZpT_[skim])->GetBinContent(j),2) +
          ((TH1D*)trig_unc_ZpT_[skim])->GetBinContent(j) +
          ((TH1D*)stat_sf_ZpT_[skim])->GetBinContent(j) +
          ((TH1D*)stat_events_ZpT_[skim])->GetBinContent(j) 
        ) /
        ((TH1D*)dN_dZpT_[skim])->GetBinContent(j)
      );
      ((TH1D*)trig_unc_ZpT_[skim])->SetBinContent(j,
        sqrt(((TH1D*)trig_unc_ZpT_[skim])->GetBinContent(j)) / ((TH1D*)dN_dZpT_[skim])->GetBinContent(j) );
      ((TH1D*)stat_sf_ZpT_[skim])->SetBinContent(j,
        sqrt(((TH1D*)stat_sf_ZpT_[skim])->GetBinContent(j)) / ((TH1D*)dN_dZpT_[skim])->GetBinContent(j) );
      ((TH1D*)stat_events_ZpT_[skim])->SetBinContent(j,
        sqrt(((TH1D*)stat_events_ZpT_[skim])->GetBinContent(j)) / ((TH1D*)dN_dZpT_[skim])->GetBinContent(j) );
      ((TH1D*)dN_dZpT_[skim])->SetBinError(j, 
        ((TH1D*)dN_dZpT_[skim])->GetBinContent(j) *  ((TH1D*)unc_ZpT_[skim])->GetBinContent(j) 
      );
      

    }
    ((TH1D*)syst_sf_ZpT_[skim])->Divide(((TH1D*)dN_dZpT_[skim]));
  }
  TClonesArray c_dN_dZpT_ ("TCanvas", n_skims);
  new(c_dN_dZpT_[0]) TCanvas("c_mc_dN_dZpT_ele_tight",   "p_{T}^{ll} distribution of Z #rightarrow ee events in MC");
  new(c_dN_dZpT_[1]) TCanvas("c_mc_dN_dZpT_mu_tight",    "p_{T}^{ll} distribution of Z #rightarrow #mu#mu events in MC");
  new(c_dN_dZpT_[2]) TCanvas("c_data_dN_dZpT_ele_tight", "p_{T}^{ll} distribution of Z #rightarrow ee events in data");
  new(c_dN_dZpT_[3]) TCanvas("c_data_dN_dZpT_mu_tight",  "p_{T}^{ll} distribution of Z #rightarrow #mu#mu events in data");
  TClonesArray c_syst_sf_ZpT_ ("TCanvas", n_skims);
  new(c_syst_sf_ZpT_[0]) TCanvas("c_mc_syst_sf_ZpT_ele_tight",   "S.F. syst. uncertainty as p_{T} of Z #rightarrow ee events in MC");
  new(c_syst_sf_ZpT_[1]) TCanvas("c_mc_syst_sf_ZpT_mu_tight",    "S.F. syst. uncertainty as p_{T} of Z #rightarrow #mu#mu events in MC");
  new(c_syst_sf_ZpT_[2]) TCanvas("c_data_syst_sf_ZpT_ele_tight", "S.F. syst. uncertainty as p_{T} of Z #rightarrow ee events in data");
  new(c_syst_sf_ZpT_[3]) TCanvas("c_data_syst_sf_ZpT_mu_tight",  "S.F. syst. uncertainty as p_{T} of Z #rightarrow #mu#mu events in data");
  TClonesArray c_stat_sf_ZpT_ ("TCanvas", n_skims);
  new(c_stat_sf_ZpT_[0]) TCanvas("c_mc_stat_sf_ZpT_ele_tight",   "S.F. stat. uncertainty as p_{T} of Z #rightarrow ee events in MC");
  new(c_stat_sf_ZpT_[1]) TCanvas("c_mc_stat_sf_ZpT_mu_tight",    "S.F. stat. uncertainty as p_{T} of Z #rightarrow #mu#mu events in MC");
  new(c_stat_sf_ZpT_[2]) TCanvas("c_data_stat_sf_ZpT_ele_tight", "S.F. stat. uncertainty as p_{T} of Z #rightarrow ee events in data");
  new(c_stat_sf_ZpT_[3]) TCanvas("c_data_stat_sf_ZpT_mu_tight",  "S.F. stat. uncertainty as p_{T} of Z #rightarrow #mu#mu events in data");
  TClonesArray c_trig_unc_ZpT_ ("TCanvas", n_skims);
  new(c_trig_unc_ZpT_[0]) TCanvas("c_mc_trig_unc_ZpT_ele_tight",   "Trigger uncertainty as p_{T} of Z #rightarrow ee events in MC");
  new(c_trig_unc_ZpT_[1]) TCanvas("c_mc_trig_unc_ZpT_mu_tight",    "Trigger uncertainty as p_{T} of Z #rightarrow #mu#mu events in MC");
  new(c_trig_unc_ZpT_[2]) TCanvas("c_data_trig_unc_ZpT_ele_tight", "Trigger uncertainty as p_{T} of Z #rightarrow ee events in data");
  new(c_trig_unc_ZpT_[3]) TCanvas("c_data_trig_unc_ZpT_mu_tight",  "Trigger uncertainty as p_{T} of Z #rightarrow #mu#mu events in data");
  TClonesArray c_unc_ZpT_ ("TCanvas", n_skims);
  new(c_unc_ZpT_[0]) TCanvas("c_mc_unc_ZpT_ele_tight",   "Yield rel. uncertainty as p_{T} of Z #rightarrow ee events in MC");
  new(c_unc_ZpT_[1]) TCanvas("c_mc_unc_ZpT_mu_tight",    "Yield rel. uncertainty as p_{T} of Z #rightarrow #mu#mu events in MC");
  new(c_unc_ZpT_[2]) TCanvas("c_data_unc_ZpT_ele_tight", "Yield rel. uncertainty as p_{T} of Z #rightarrow ee events in data");
  new(c_unc_ZpT_[3]) TCanvas("c_data_unc_ZpT_mu_tight",  "Yield rel. uncertainty as p_{T} of Z #rightarrow #mu#mu events in data");
  
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;
  Int_t mit_red  = 1861; 
  Int_t mit_gray = 1862; 
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);
  
  char histo_file_name[256];
  for(int skim=0; skim<n_skims; skim++) {
    
    //styling
    ((TCanvas*)c_dN_dZpT_[skim])->cd();
    ((TCanvas*)c_dN_dZpT_[skim])->SetLogy();
    ((TCanvas*)c_dN_dZpT_[skim])->SetLogx();
    ( (TH1D*) dN_dZpT_[skim]    )->Draw();
    sprintf(histo_file_name, "%s%s.png", plots_dir.c_str(), ((TH1D*)dN_dZpT_[skim])->GetName());
    ((TCanvas*)c_dN_dZpT_[skim])->Print(histo_file_name);
    if(skim<=1) {
      ((TCanvas*)c_syst_sf_ZpT_[skim])->cd();
      ((TCanvas*)c_syst_sf_ZpT_[skim])->SetLogx();
      ( (TH1D*) syst_sf_ZpT_[skim]    )->Draw("LP");
      sprintf(histo_file_name, "%s%s.png", plots_dir.c_str(), ((TH1D*)syst_sf_ZpT_[skim])->GetName());
      ((TCanvas*)c_syst_sf_ZpT_[skim])->Print(histo_file_name);
      ((TCanvas*)c_stat_sf_ZpT_[skim])->cd();
      ((TCanvas*)c_stat_sf_ZpT_[skim])->SetLogx();
      ( (TH1D*) stat_sf_ZpT_[skim]    )->Draw("LP");
      sprintf(histo_file_name, "%s%s.png", plots_dir.c_str(), ((TH1D*)stat_sf_ZpT_[skim])->GetName());
      ((TCanvas*)c_stat_sf_ZpT_[skim])->Print(histo_file_name);
      ((TCanvas*)c_trig_unc_ZpT_[skim])->SetLogx();
      ((TCanvas*)c_trig_unc_ZpT_[skim])->cd();
      ( (TH1D*) trig_unc_ZpT_[skim]    )->Draw("LP");
      sprintf(histo_file_name, "%s%s.png", plots_dir.c_str(), ((TH1D*)trig_unc_ZpT_[skim])->GetName());
      ((TCanvas*)c_trig_unc_ZpT_[skim])->Print(histo_file_name);
    }
    ((TCanvas*)c_unc_ZpT_[skim])->cd();
    //( (TH1D*) unc_ZpT_[skim]    )->Divide( (TH1D*) dN_dZpT_[skim] );
    ( (TH1D*) unc_ZpT_[skim]    )->Draw("LP");
    sprintf(histo_file_name, "%s%s.png", plots_dir.c_str(), ((TH1D*)unc_ZpT_[skim])->GetName());
    ((TCanvas*)c_unc_ZpT_[skim])->Print(histo_file_name);
  }
  
  TH1D *flavor_ratio_mc   = new TH1D("flavor_ratio_mc",   "", nbins, Z_pT_bins);
  TH1D *flavor_ratio_data = new TH1D("flavor_ratio_data", "", nbins, Z_pT_bins);
  (*flavor_ratio_mc)   = (* ((TH1D*)dN_dZpT_[1])) / (* ((TH1D*)dN_dZpT_[0]));
  (*flavor_ratio_data) = (* ((TH1D*)dN_dZpT_[3])) / (* ((TH1D*)dN_dZpT_[2]));
  flavor_ratio_mc->SetTitle("Ratio of Z #rightarrow #mu#mu / Z #rightarrow ee events");
  //flavor_ratio_mc->SetFillStyle(3002);
  //flavor_ratio_mc->SetFillColor(mit_red);
  flavor_ratio_mc->SetLineColor(mit_red);
  flavor_ratio_mc->SetMarkerColor(mit_red);
  flavor_ratio_mc->SetMarkerStyle(20);
  flavor_ratio_mc->SetMaximum(3);
  flavor_ratio_mc->SetMinimum(1);
  flavor_ratio_data->SetTitle("Ratio of Z #rightarrow #mu#mu / Z #rightarrow ee events");
  //flavor_ratio_data->SetFillStyle(3004);
  //flavor_ratio_data->SetFillColor(mit_gray);
  flavor_ratio_data->SetLineColor(mit_gray);
  flavor_ratio_data->SetMarkerColor(1);
  flavor_ratio_data->SetMarkerStyle(21);
  //TCanvas *c_flavor_ratio=new TCanvas("c_flavor_ratio", "c_flavor_ratio"); 
  //c_flavor_ratio->SetLogx();
  //c_flavor_ratio->SetGrid(0,1);
  //flavor_ratio_mc->SetMinimum(1);
  //flavor_ratio_mc->SetMaximum(3);
  //flavor_ratio_mc->Draw("E1 P0 L");
  //flavor_ratio_mc->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  //flavor_ratio_data->Draw("E1 P0 L SAME");
  //TLegend *l_flavor_ratio = new TLegend(.7,.7,.85,.85);
  //l_flavor_ratio->AddEntry(flavor_ratio_mc, "MC", "lp");
  //l_flavor_ratio->AddEntry(flavor_ratio_data, "Data (2.1 fb^{-1})", "lp");
  //l_flavor_ratio->SetFillColor(0);
  //l_flavor_ratio->Draw("SAME");
  //c_flavor_ratio->Print((plots_dir+"ZpT_flavor_ratio.png").c_str());
   
  TCanvas *c_flavor_ratio = ratio_plot(
    "flavor_ratio",
    "Ratio of Z #rightarrow #mu#mu / Z #rightarrow ee events",
    "Z p_{T} [GeV]",
    flavor_ratio_data,
    flavor_ratio_mc,
    true,
    false,
    false
  );
 
  c_flavor_ratio->Print((plots_dir+"ZpT_flavor_ratio.png").c_str());
  
  printf("Done! \n");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void electron_closure(
  string plots_dir,
  string root_dir,
  bool verbose=false
) {
  // Pad directories with a slash at the end if it's not there
  if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
  if( root_dir[plots_dir.size()-1]  != '/' )  root_dir  = root_dir + "/";
  gStyle->SetOptStat(0); 
  int n_skims=4; 
  int nbins=8;
  Float_t Z_pT_bins[] = {50., 75., 100., 125., 150., 175., 200., 300., 400.};
  int max_pT = 400;
  TFile *f_sf_ele  = TFile::Open((root_dir+"scalefactors_ele_74x.root").c_str(),"READ");
  TFile *f_ele_triggers = TFile::Open("/home/dhsu/TagAndProbe/Data_Triggers/Ele27_eta2p1_WPLoose_Gsf/eff.root", "READ");
  
  TFile *f_ele_mc = TFile::Open("~/leptonScaleFactors/root/DYJetsToLL_74x_aMCatNLO_genMatching_BaselineToTight_electronTnP.root","READ");
  TFile *f_ele_data = TFile::Open("~/leptonScaleFactors/root/SingleElectron_BaselineToTight_electronTnP.root","READ");

  TH2D *syst_ele_sf_tight    = (TH2D*) f_sf_ele ->Get("scalefactors_Tight_ele_syst_error_combined");
  TH2D *stat_ele_sf_tight_lo = (TH2D*) f_sf_ele ->Get("scalefactors_Tight_ele_stat_error_lo");
  TH2D *stat_ele_sf_tight_hi = (TH2D*) f_sf_ele ->Get("scalefactors_Tight_ele_stat_error_hi");

  TH2D *trig_eff_ele = (TH2D*) f_ele_triggers->Get("hEffEtaPt");
  TH2D *unc_trig_eff_ele_hi = (TH2D*) f_ele_triggers->Get("hErrhEtaPt");
  TH2D *unc_trig_eff_ele_lo = (TH2D*) f_ele_triggers->Get("hErrlEtaPt");
  TH2D *sf_ele_tight = (TH2D*) f_sf_ele->Get("scalefactors_Tight_ele");
  
  // read from tnp skim
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
  
  TTree *mc_tree = (TTree*) f_ele_mc->Get("Events");
  mc_tree->SetBranchAddress("runNum",   &runNum   );  
  mc_tree->SetBranchAddress("lumiSec",  &lumiSec  );  
  mc_tree->SetBranchAddress("evtNum",   &evtNum   );  
  mc_tree->SetBranchAddress("npv",      &npv      );  
  mc_tree->SetBranchAddress("pass",     &pass     );  
  mc_tree->SetBranchAddress("npu",      &npu      );  
  mc_tree->SetBranchAddress("scale1fb", &scale1fb );
  mc_tree->SetBranchAddress("mass",     &mass     );  
  mc_tree->SetBranchAddress("qtag",     &qtag     );  
  mc_tree->SetBranchAddress("qprobe",   &qprobe   );  
  mc_tree->SetBranchAddress("tag",      &p4_tag   );  
  mc_tree->SetBranchAddress("probe",    &p4_probe );     
  TTree *data_tree = (TTree*) f_ele_data->Get("Events");
  data_tree->SetBranchAddress("runNum",   &runNum   );  
  data_tree->SetBranchAddress("lumiSec",  &lumiSec  );  
  data_tree->SetBranchAddress("evtNum",   &evtNum   );  
  data_tree->SetBranchAddress("npv",      &npv      );  
  data_tree->SetBranchAddress("pass",     &pass     );  
  data_tree->SetBranchAddress("npu",      &npu      );  
  data_tree->SetBranchAddress("scale1fb", &scale1fb );
  data_tree->SetBranchAddress("mass",     &mass     );  
  data_tree->SetBranchAddress("qtag",     &qtag     );  
  data_tree->SetBranchAddress("qprobe",   &qprobe   );  
  data_tree->SetBranchAddress("tag",      &p4_tag   );  
  data_tree->SetBranchAddress("probe",    &p4_probe );     
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);
  
  TH1D *probe_pt_mc_EBEB = new TH1D("probe_pt_mc_EBEB", "", 95,10,200);
  probe_pt_mc_EBEB->Sumw2(kTRUE);
  probe_pt_mc_EBEB->SetDefaultSumw2(kTRUE);
  TH1D *probe_pt_mc_EEEE = new TH1D("probe_pt_mc_EEEE", "", 95,10,200);
  TH1D *tag_pt_mc_EBEB = new TH1D("tag_pt_mc_EBEB", "", 85,30,200);
  TH1D *tag_pt_mc_EEEE = new TH1D("tag_pt_mc_EEEE", "", 85,30,200);
  TH1D *probe_pt_data_EBEB = new TH1D("probe_pt_data_EBEB", "", 95,10,200);
  TH1D *probe_pt_data_EEEE = new TH1D("probe_pt_data_EEEE", "", 95,10,200);
  TH1D *tag_pt_data_EBEB = new TH1D("tag_pt_data_EBEB", "", 85,30,200);
  TH1D *tag_pt_data_EEEE = new TH1D("tag_pt_data_EEEE", "", 85,30,200);

  TH1D *tag_eta_mc   = new TH1D("tag_eta_mc", "", 42, -2.1, 2.1);
  TH1D *tag_eta_data   = new TH1D("tag_eta_data", "", 42, -2.1, 2.1);
  TH1D *probe_eta_mc = new TH1D("probe_eta_mc", "", 50, -2.5, 2.5);
  TH1D *probe_eta_data = new TH1D("probe_eta_data", "", 50, -2.5, 2.5);

  TH1D *dilepton_pt_mc_EBEB = new TH1D("dilepton_pt_mc_EBEB", "", 100,0,200);
  TH1D *dilepton_pt_mc_EBEE = new TH1D("dilepton_pt_mc_EBEE", "", 100,0,200);
  TH1D *dilepton_pt_mc_EEEE = new TH1D("dilepton_pt_mc_EEEE", "", 100,0,200);
  TH1D *dilepton_pt_data_EBEB = new TH1D("dilepton_pt_data_EBEB", "", 100,0,200);
  TH1D *dilepton_pt_data_EBEE = new TH1D("dilepton_pt_data_EBEE", "", 100,0,200);
  TH1D *dilepton_pt_data_EEEE = new TH1D("dilepton_pt_data_EEEE", "", 100,0,200);

  TH1D *dilepton_M_mc_EBEB = new TH1D("dilepton_M_mc_EBEB", "", 60,60,120);
  TH1D *dilepton_M_mc_EBEE = new TH1D("dilepton_M_mc_EBEE", "", 60,60,120);
  TH1D *dilepton_M_mc_EEEE = new TH1D("dilepton_M_mc_EEEE", "", 60,60,120);
  TH1D *dilepton_M_data_EBEB = new TH1D("dilepton_M_data_EBEB", "", 60,60,120);
  TH1D *dilepton_M_data_EBEE = new TH1D("dilepton_M_data_EBEE", "", 60,60,120);
  TH1D *dilepton_M_data_EEEE = new TH1D("dilepton_M_data_EEEE", "", 60,60,120);
  Long64_t nentries;
  double mc_reweighting = 2151*6025.2/19312436.;
  printf("mc reweighting factor is %f\n", mc_reweighting);
  // dont handle mc xs error for now
  nentries = mc_tree->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    mc_tree->GetEntry(i);
    if(!(
      pass==1 &&
      qtag + qprobe == 0 && 
      TMath::Abs(mass-90.) <= 30
    )) continue;
    TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
    double Z_pT = systemP4.Pt();
    double trig_eff_tag=1, trig_eff_probe=1, sf_tag=1, sf_probe=1;
    double stat_sf_tag=0, stat_sf_probe=0, syst_sf_tag=0, syst_sf_probe=0, unc_trig_tag=0, unc_trig_probe;
    double tag_eta = p4_tag->Eta();
    double tag_pT = p4_tag->Pt();
    double probe_eta = p4_probe->Eta();
    double probe_pT = p4_probe->Pt();
    if(tag_pT < 30) continue;
    if(tag_eta >= 2.1)    continue;
    if(tag_eta <= -2.1)   continue;
    if(probe_eta >= 2.5)  continue;
    if(probe_eta <= -2.5) continue;
    sf_tag = sf_ele_tight->GetBinContent( sf_ele_tight->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT )));
    sf_probe = sf_ele_tight->GetBinContent( sf_ele_tight->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT )));
    
    stat_sf_tag   = TMath::Max(
      stat_ele_sf_tight_lo->GetBinContent(stat_ele_sf_tight_lo->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT))),
      stat_ele_sf_tight_hi->GetBinContent(stat_ele_sf_tight_hi->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT)))
    ); 
    stat_sf_probe = TMath::Max( 
      stat_ele_sf_tight_lo->GetBinContent(stat_ele_sf_tight_lo->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT))),
      stat_ele_sf_tight_hi->GetBinContent(stat_ele_sf_tight_hi->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT)))
    );
    syst_sf_tag   = syst_ele_sf_tight->GetBinContent(syst_ele_sf_tight->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT)));
    syst_sf_probe = syst_ele_sf_tight->GetBinContent(syst_ele_sf_tight->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT)));
    
    trig_eff_tag   = trig_eff_ele->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT) ) );
    if(probe_pT < 30) trig_eff_probe = 0;
    else              trig_eff_probe = trig_eff_ele->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT) ) );

    unc_trig_tag   = TMath::Max(
      unc_trig_eff_ele_lo->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(tag_eta),  tag_pT ) ),
      unc_trig_eff_ele_hi->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(tag_eta),  tag_pT ) )
    );
    if(probe_pT < 30) unc_trig_probe = 0;
    else              unc_trig_probe = TMath::Max(
      unc_trig_eff_ele_hi->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT) )),
      unc_trig_eff_ele_lo->GetBinContent( trig_eff_ele->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT) ))
    );
    //double weight = scale1fb*mc_reweighting;
    double weight = sf_tag * sf_probe * (1 - (1 - trig_eff_tag)*(1 - trig_eff_probe)) * scale1fb * mc_reweighting;

    // EB EB
    if(TMath::Abs(tag_eta) < 1.4442 && TMath::Abs(probe_eta) < 1.4442) {
      probe_pt_mc_EBEB    ->Fill(probe_pT,  weight);
      probe_eta_mc        ->Fill(probe_eta, weight);
      tag_pt_mc_EBEB      ->Fill(tag_pT,    weight);
      tag_eta_mc          ->Fill(tag_eta,   weight);
      dilepton_M_mc_EBEB  ->Fill(mass,      weight);
      dilepton_pt_mc_EBEB ->Fill(Z_pT,      weight);
    }
    // EB EE
    else if(
      (TMath::Abs(tag_eta) < 1.4442 && TMath::Abs(probe_eta) > 1.566) ||
      (TMath::Abs(tag_eta) > 1.566  && TMath::Abs(probe_eta) < 1.4442)
    ) {
      probe_eta_mc        ->Fill(probe_eta, weight);
      tag_eta_mc          ->Fill(tag_eta,   weight);
      dilepton_M_mc_EBEE  ->Fill(mass,      weight);
      dilepton_pt_mc_EBEE ->Fill(Z_pT,      weight);
      
    }
    // EE EE
    else if(TMath::Abs(tag_eta) > 1.566  && TMath::Abs(probe_eta) > 1.566) {
      probe_pt_mc_EEEE    ->Fill(probe_pT,  weight);
      probe_eta_mc        ->Fill(probe_eta, weight);
      tag_pt_mc_EEEE      ->Fill(tag_pT,    weight);
      tag_eta_mc          ->Fill(tag_eta,   weight);
      dilepton_M_mc_EEEE  ->Fill(mass,      weight);
      dilepton_pt_mc_EEEE ->Fill(Z_pT,      weight);

    }

  }
  printf("done reading MC skim\n");
  nentries = data_tree->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    data_tree->GetEntry(i);
    if(!(
      pass==1 &&
      qtag + qprobe == 0 &&
      TMath::Abs(mass-90.) <= 30
    )) continue;
    TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
    double Z_pT = systemP4.Pt();
    double tag_eta = p4_tag->Eta();
    double tag_pT = p4_tag->Pt();
    double probe_eta = p4_probe->Eta();
    double probe_pT = p4_probe->Pt();
    if(tag_pT < 30) continue;
    if(tag_eta >= 2.1)    continue;
    if(tag_eta <= -2.1)   continue;
    if(probe_eta >= 2.5)  continue;
    if(probe_eta <= -2.5) continue;

    double weight=1;
    // EB EB
    if(TMath::Abs(tag_eta) < 1.4442 && TMath::Abs(probe_eta) < 1.4442) {
      probe_pt_data_EBEB    ->Fill(probe_pT,  weight);
      probe_eta_data        ->Fill(probe_eta, weight);
      tag_pt_data_EBEB      ->Fill(tag_pT,    weight);
      tag_eta_data          ->Fill(tag_eta,   weight);
      dilepton_M_data_EBEB  ->Fill(mass,      weight);
      dilepton_pt_data_EBEB ->Fill(Z_pT,      weight);
    }
    // EB EE
    else if(
      (TMath::Abs(tag_eta) < 1.4442 && TMath::Abs(probe_eta) > 1.566) ||
      (TMath::Abs(tag_eta) > 1.566  && TMath::Abs(probe_eta) < 1.4442)
    ) {
      probe_eta_data        ->Fill(probe_eta, weight);
      tag_eta_data          ->Fill(tag_eta,   weight);
      dilepton_M_data_EBEE  ->Fill(mass,      weight);
      dilepton_pt_data_EBEE ->Fill(Z_pT,      weight);
      
    }
    // EE EE
    else if(TMath::Abs(tag_eta) > 1.566  && TMath::Abs(probe_eta) > 1.566) {
      probe_pt_data_EEEE    ->Fill(probe_pT,  weight);
      probe_eta_data        ->Fill(probe_eta, weight);
      tag_pt_data_EEEE      ->Fill(tag_pT,    weight);
      tag_eta_data          ->Fill(tag_eta,   weight);
      dilepton_M_data_EEEE  ->Fill(mass,      weight);
      dilepton_pt_data_EEEE ->Fill(Z_pT,      weight);

    }

  }
  printf("done reading data skim\n");
  //mc_tree  ->Draw("probe.Pt() >> probe_pt_mc_EBEB", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() < 1.4442 && tag.Eta() < 1.4442)");
  //mc_tree  ->Draw("probe.Pt() >> probe_pt_mc_EEEE", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() > 1.566 && tag.Eta() > 1.566)");
  //mc_tree  ->Draw("tag.Pt() >> tag_pt_mc_EBEB", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() < 1.4442 && tag.Eta() < 1.4442)");
  //mc_tree  ->Draw("tag.Pt() >> tag_pt_mc_EEEE", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() > 1.566 && tag.Eta() > 1.566)");
  //data_tree->Draw("probe.Pt() >> probe_pt_data_EBEB", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() < 1.4442 && tag.Eta() < 1.4442)");
  //data_tree->Draw("probe.Pt() >> probe_pt_data_EEEE", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() > 1.566 && tag.Eta() > 1.566)");
  //data_tree->Draw("tag.Pt()   >> tag_pt_data_EBEB",   "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() < 1.4442 && tag.Eta() < 1.4442)");
  //data_tree->Draw("tag.Pt()   >> tag_pt_data_EEEE",   "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() > 1.566 && tag.Eta() > 1.566)");
  //data_tree->Draw("probe.Eta() >> probe_eta_data", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30)");
  //data_tree->Draw("tag.Eta() >> tag_eta_data", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30)");
  //
  //data_tree->Draw("(probe+tag).Pt() >> dilepton_pt_data_EBEB", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() < 1.4442 && tag.Eta() < 1.4442)");
  //data_tree->Draw("(probe+tag).Pt() >> dilepton_pt_data_EBEE", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && ((probe.Eta() < 1.4442 && tag.Eta() > 1.566) || (probe.Eta() > 1.566 && tag.Eta() < 1.4442)))");
  //data_tree->Draw("(probe+tag).Pt() >> dilepton_pt_data_EEEE", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() > 1.566 && tag.Eta() > 1.566)");
  //
  //data_tree->Draw("(probe+tag).M() >> dilepton_M_data_EBEB", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() < 1.4442 && tag.Eta() < 1.4442)");
  //data_tree->Draw("(probe+tag).M() >> dilepton_M_data_EBEE", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && ((probe.Eta() < 1.4442 && tag.Eta() > 1.566) || (probe.Eta() > 1.566 && tag.Eta() < 1.4442)))");
  //data_tree->Draw("(probe+tag).M() >> dilepton_M_data_EEEE", "scale1fb*(pass==1 && (qtag + qprobe == 0) && tag.Pt() > 30 && probe.Eta() > 1.566 && tag.Eta() > 1.566)");
  
  TCanvas *Zee_probe_pt_EBEB    = ratio_plot("probe_pt_EBEB", "Probe p_{T} (EB-EB)", "p_{T} [GeV]", probe_pt_data_EBEB, probe_pt_mc_EBEB);
  TCanvas *Zee_probe_pt_EEEE    = ratio_plot("probe_pt_EEEE", "Probe p_{T} (EE-EE)", "p_{T} [GeV]", probe_pt_data_EEEE, probe_pt_mc_EEEE);
  TCanvas *Zee_tag_pt_EBEB      = ratio_plot("tag_pt_EBEB",   "Tag p_{T} (EB-EB)", "p_{T} [GeV]", tag_pt_data_EBEB, tag_pt_mc_EBEB);
  TCanvas *Zee_tag_pt_EEEE      = ratio_plot("tag_pt_EEEE",   "Tag p_{T} (EE-EE)", "p_{T} [GeV]", tag_pt_data_EEEE, tag_pt_mc_EEEE);
  TCanvas *Zee_tag_eta          = ratio_plot("tag_eta",   "Tag #eta", "#eta", tag_eta_data, tag_eta_mc);
  TCanvas *Zee_probe_eta        = ratio_plot("probe_eta", "Probe #eta", "#eta", probe_eta_data, probe_eta_mc);
  TCanvas *Zee_dilepton_pt_EBEB = ratio_plot("dilepton_pt_EBEB", "p_{T}^{ll} (EB-EB)", "p_{T} [GeV]", dilepton_pt_data_EBEB, dilepton_pt_mc_EBEB);
  TCanvas *Zee_dilepton_pt_EBEE = ratio_plot("dilepton_pt_EBEE", "p_{T}^{ll} (EB-EE)", "p_{T} [GeV]", dilepton_pt_data_EBEE, dilepton_pt_mc_EBEE);
  TCanvas *Zee_dilepton_pt_EEEE = ratio_plot("dilepton_pt_EEEE", "p_{T}^{ll} (EE-EE)", "p_{T} [GeV]", dilepton_pt_data_EEEE, dilepton_pt_mc_EEEE);
  TCanvas *Zee_dilepton_M_EBEB  = ratio_plot("dilepton_M_EBEB",  "m^{ll} (EB-EB)", "invariant mass [GeV]", dilepton_M_data_EBEB, dilepton_M_mc_EBEB);
  TCanvas *Zee_dilepton_M_EBEE  = ratio_plot("dilepton_M_EBEE",  "m^{ll} (EB-EE)", "invariant mass [GeV]", dilepton_M_data_EBEE, dilepton_M_mc_EBEE);
  TCanvas *Zee_dilepton_M_EEEE  = ratio_plot("dilepton_M_EEEE",  "m^{ll} (EE-EE)", "invariant mass [GeV]", dilepton_M_data_EEEE, dilepton_M_mc_EEEE);

  Zee_probe_pt_EBEB    ->Print( (plots_dir + string(Zee_probe_pt_EBEB    ->GetName()) + ".png").c_str());
  Zee_probe_pt_EEEE    ->Print( (plots_dir + string(Zee_probe_pt_EEEE    ->GetName()) + ".png").c_str());
  Zee_tag_pt_EBEB      ->Print( (plots_dir + string(Zee_tag_pt_EBEB      ->GetName()) + ".png").c_str());
  Zee_tag_pt_EEEE      ->Print( (plots_dir + string(Zee_tag_pt_EEEE      ->GetName()) + ".png").c_str());
  Zee_tag_eta          ->Print( (plots_dir + string(Zee_tag_eta          ->GetName()) + ".png").c_str());
  Zee_probe_eta        ->Print( (plots_dir + string(Zee_probe_eta        ->GetName()) + ".png").c_str());
  Zee_dilepton_pt_EBEB ->Print( (plots_dir + string(Zee_dilepton_pt_EBEB ->GetName()) + ".png").c_str());
  Zee_dilepton_pt_EBEE ->Print( (plots_dir + string(Zee_dilepton_pt_EBEE ->GetName()) + ".png").c_str());
  Zee_dilepton_pt_EEEE ->Print( (plots_dir + string(Zee_dilepton_pt_EEEE ->GetName()) + ".png").c_str());
  Zee_dilepton_M_EBEB  ->Print( (plots_dir + string(Zee_dilepton_M_EBEB  ->GetName()) + ".png").c_str());
  Zee_dilepton_M_EBEE  ->Print( (plots_dir + string(Zee_dilepton_M_EBEE  ->GetName()) + ".png").c_str());
  Zee_dilepton_M_EEEE  ->Print( (plots_dir + string(Zee_dilepton_M_EEEE  ->GetName()) + ".png").c_str());


   
  printf("Done! \n");
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void muon_closure(
  string plots_dir,
  string root_dir,
  bool verbose=false
) {
  // Pad directories with a slash at the end if it's not there
  if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
  if( root_dir[plots_dir.size()-1]  != '/' )  root_dir  = root_dir + "/";
  gStyle->SetOptStat(0); 
  int n_skims=4; 
  int nbins=8;
  Float_t Z_pT_bins[] = {50., 75., 100., 125., 150., 175., 200., 300., 400.};
  int max_pT = 400;
  TFile *f_sf_mu  = TFile::Open((root_dir+"scalefactors_mu_74x.root").c_str(),"READ");
  TFile *f_mu_triggers = TFile::Open("/home/dhsu/TagAndProbe/Data_Triggers/IsoMu27/eff.root", "READ");
  
  TFile *f_mu_mc = TFile::Open("~/leptonScaleFactors/root/DYJetsToLL_74x_aMCatNLO_genMatching_BaselineToTight_muonTnP.root","READ");
  TFile *f_mu_data = TFile::Open("~/leptonScaleFactors/root/SingleMuon_BaselineToTight_muonTnP.root","READ");

  TH2D *syst_mu_sf_tight    = (TH2D*) f_sf_mu ->Get("scalefactors_Tight_mu_syst_error_combined");
  TH2D *stat_mu_sf_tight_lo = (TH2D*) f_sf_mu ->Get("scalefactors_Tight_mu_stat_error_lo");
  TH2D *stat_mu_sf_tight_hi = (TH2D*) f_sf_mu ->Get("scalefactors_Tight_mu_stat_error_hi");

  TH2D *trig_eff_mu = (TH2D*) f_mu_triggers->Get("hEffEtaPt");
  TH2D *unc_trig_eff_mu_hi = (TH2D*) f_mu_triggers->Get("hErrhEtaPt");
  TH2D *unc_trig_eff_mu_lo = (TH2D*) f_mu_triggers->Get("hErrlEtaPt");
  TH2D *sf_mu_tight = (TH2D*) f_sf_mu->Get("scalefactors_Tight_mu");
  
  // read from tnp skim
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
  
  TTree *mc_tree = (TTree*) f_mu_mc->Get("Events");
  mc_tree->SetBranchAddress("runNum",   &runNum   );  
  mc_tree->SetBranchAddress("lumiSec",  &lumiSec  );  
  mc_tree->SetBranchAddress("evtNum",   &evtNum   );  
  mc_tree->SetBranchAddress("npv",      &npv      );  
  mc_tree->SetBranchAddress("pass",     &pass     );  
  mc_tree->SetBranchAddress("npu",      &npu      );  
  mc_tree->SetBranchAddress("scale1fb", &scale1fb );
  mc_tree->SetBranchAddress("mass",     &mass     );  
  mc_tree->SetBranchAddress("qtag",     &qtag     );  
  mc_tree->SetBranchAddress("qprobe",   &qprobe   );  
  mc_tree->SetBranchAddress("tag",      &p4_tag   );  
  mc_tree->SetBranchAddress("probe",    &p4_probe );     
  TTree *data_tree = (TTree*) f_mu_data->Get("Events");
  data_tree->SetBranchAddress("runNum",   &runNum   );  
  data_tree->SetBranchAddress("lumiSec",  &lumiSec  );  
  data_tree->SetBranchAddress("evtNum",   &evtNum   );  
  data_tree->SetBranchAddress("npv",      &npv      );  
  data_tree->SetBranchAddress("pass",     &pass     );  
  data_tree->SetBranchAddress("npu",      &npu      );  
  data_tree->SetBranchAddress("scale1fb", &scale1fb );
  data_tree->SetBranchAddress("mass",     &mass     );  
  data_tree->SetBranchAddress("qtag",     &qtag     );  
  data_tree->SetBranchAddress("qprobe",   &qprobe   );  
  data_tree->SetBranchAddress("tag",      &p4_tag   );  
  data_tree->SetBranchAddress("probe",    &p4_probe );     
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette();
  TPaletteAxis *palette_axis;
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);
  
  TH1D *probe_pt_mc      = new TH1D("probe_pt_mc     ", "", 95,10,200);
  probe_pt_mc     ->Sumw2(kTRUE);
  probe_pt_mc     ->SetDefaultSumw2(kTRUE);
  TH1D *tag_pt_mc      = new TH1D("tag_pt_mc     ", "", 85,30,200);
  TH1D *probe_pt_data      = new TH1D("probe_pt_data     ", "", 95,10,200);
  TH1D *tag_pt_data      = new TH1D("tag_pt_data     ", "", 85,30,200);

  TH1D *tag_eta_mc   = new TH1D("tag_eta_mc", "", 42, -2.1, 2.1);
  TH1D *tag_eta_data   = new TH1D("tag_eta_data", "", 42, -2.1, 2.1);
  TH1D *probe_eta_mc = new TH1D("probe_eta_mc", "", 48, -2.4, 2.4);
  TH1D *probe_eta_data = new TH1D("probe_eta_data", "", 48, -2.4, 2.4);

  TH1D *dilepton_pt_mc      = new TH1D("dilepton_pt_mc     ", "", 100,0,200);
  TH1D *dilepton_pt_data      = new TH1D("dilepton_pt_data     ", "", 100,0,200);

  TH1D *dilepton_M_mc      = new TH1D("dilepton_M_mc     ", "", 60,60,120);
  TH1D *dilepton_M_data      = new TH1D("dilepton_M_data     ", "", 60,60,120);
  Long64_t nentries;
  //double mc_reweighting = 2150.*6025.2/19312436.;
  double mc_reweighting = (double) data_tree->GetEntries() / (double)mc_tree->GetEntries();
  printf("mc reweighting factor is %f\n", mc_reweighting);
  // dont handle mc xs error for now
  nentries = mc_tree->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    mc_tree->GetEntry(i);
    if(!(
      pass==1 &&
      qtag + qprobe == 0 && 
      TMath::Abs(mass-90.) <= 30
    )) continue;
    TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
    double Z_pT = systemP4.Pt();
    double trig_eff_tag=1, trig_eff_probe=1, sf_tag=1, sf_probe=1;
    double stat_sf_tag=0, stat_sf_probe=0, syst_sf_tag=0, syst_sf_probe=0, unc_trig_tag=0, unc_trig_probe;
    double tag_eta = p4_tag->Eta();
    double tag_pT = p4_tag->Pt();
    double probe_eta = p4_probe->Eta();
    double probe_pT = p4_probe->Pt();
    if(tag_pT < 30) continue;
    if(tag_eta >= 2.1)    continue;
    if(tag_eta <= -2.1)   continue;
    if(probe_eta >= 2.4)  continue;
    if(probe_eta <= -2.4) continue;
    sf_tag = sf_mu_tight->GetBinContent( sf_mu_tight->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT )));
    sf_probe = sf_mu_tight->GetBinContent( sf_mu_tight->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT )));
    
    stat_sf_tag   = TMath::Max(
      stat_mu_sf_tight_lo->GetBinContent(stat_mu_sf_tight_lo->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT))),
      stat_mu_sf_tight_hi->GetBinContent(stat_mu_sf_tight_hi->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT)))
    ); 
    stat_sf_probe = TMath::Max( 
      stat_mu_sf_tight_lo->GetBinContent(stat_mu_sf_tight_lo->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT))),
      stat_mu_sf_tight_hi->GetBinContent(stat_mu_sf_tight_hi->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT)))
    );
    syst_sf_tag   = syst_mu_sf_tight->GetBinContent(syst_mu_sf_tight->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT)));
    syst_sf_probe = syst_mu_sf_tight->GetBinContent(syst_mu_sf_tight->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT)));
    
    trig_eff_tag   = trig_eff_mu->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(tag_eta), TMath::Min(199.9, tag_pT) ) );
    if(probe_pT < 30) trig_eff_probe = 0;
    else              trig_eff_probe = trig_eff_mu->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT) ) );

    unc_trig_tag   = TMath::Max(
      unc_trig_eff_mu_lo->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(tag_eta),  tag_pT ) ),
      unc_trig_eff_mu_hi->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(tag_eta),  tag_pT ) )
    );
    if(probe_pT < 30) unc_trig_probe = 0;
    else              unc_trig_probe = TMath::Max(
      unc_trig_eff_mu_hi->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT) )),
      unc_trig_eff_mu_lo->GetBinContent( trig_eff_mu->FindBin( TMath::Abs(probe_eta), TMath::Min(199.9, probe_pT) ))
    );
    //double weight = scale1fb*mc_reweighting;
    double weight = sf_tag * sf_probe * (1 - (1 - trig_eff_tag)*(1 - trig_eff_probe)) * scale1fb * mc_reweighting;

    probe_pt_mc         ->Fill(probe_pT,  weight);
    probe_eta_mc        ->Fill(probe_eta, weight);
    tag_pt_mc           ->Fill(tag_pT,    weight);
    tag_eta_mc          ->Fill(tag_eta,   weight);
    dilepton_M_mc       ->Fill(mass,      weight);
    dilepton_pt_mc      ->Fill(Z_pT,      weight);

  }
  printf("done reading MC skim\n");
  nentries = data_tree->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    data_tree->GetEntry(i);
    if(!(
      pass==1 &&
      qtag + qprobe == 0 &&
      TMath::Abs(mass-90.) <= 30
    )) continue;
    TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
    double Z_pT = systemP4.Pt();
    double tag_eta = p4_tag->Eta();
    double tag_pT = p4_tag->Pt();
    double probe_eta = p4_probe->Eta();
    double probe_pT = p4_probe->Pt();
    if(tag_pT < 30) continue;
    if(tag_eta >= 2.1)    continue;
    if(tag_eta <= -2.1)   continue;
    if(probe_eta >= 2.4)  continue;
    if(probe_eta <= -2.4) continue;

    double weight=1;
    probe_pt_data         ->Fill(probe_pT,  weight);
    probe_eta_data        ->Fill(probe_eta, weight);
    tag_pt_data           ->Fill(tag_pT,    weight);
    tag_eta_data          ->Fill(tag_eta,   weight);
    dilepton_M_data       ->Fill(mass,      weight);
    dilepton_pt_data      ->Fill(Z_pT,      weight);
  }
  printf("done reading data skim\n");
  
  TCanvas *Zmm_probe_pt         = ratio_plot("probe_pt", "Probe p_{T}", "p_{T} [GeV]", probe_pt_data     , probe_pt_mc     );
  TCanvas *Zmm_tag_pt           = ratio_plot("tag_pt",   "Tag p_{T}", "p_{T} [GeV]", tag_pt_data     , tag_pt_mc     );
  TCanvas *Zmm_tag_eta          = ratio_plot("tag_eta",   "Tag #eta", "#eta", tag_eta_data, tag_eta_mc);
  TCanvas *Zmm_probe_eta        = ratio_plot("probe_eta", "Probe #eta", "#eta", probe_eta_data, probe_eta_mc);
  TCanvas *Zmm_dilepton_pt      = ratio_plot("dilepton_pt", "p_{T}^{ll}", "p_{T} [GeV]", dilepton_pt_data     , dilepton_pt_mc     );
  TCanvas *Zmm_dilepton_M       = ratio_plot("dilepton_M",  "m^{ll}", "invariant mass [GeV]", dilepton_M_data     , dilepton_M_mc     );

  Zmm_probe_pt         ->Print( (plots_dir + "Zmm_probe_pt.png").c_str());
  Zmm_tag_pt           ->Print( (plots_dir + "Zmm_tag_pt.png").c_str());
  Zmm_tag_eta          ->Print( (plots_dir + "Zmm_tag_eta.png").c_str());
  Zmm_probe_eta        ->Print( (plots_dir + "Zmm_probe_eta.png").c_str());
  Zmm_dilepton_pt      ->Print( (plots_dir + "Zmm_dilepton_pt.png").c_str());
  Zmm_dilepton_M       ->Print( (plots_dir + "Zmm_dilepton_M.png").c_str());
   
  printf("Done! \n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void egamma_comparison_textfile() {
  TFile *dgh_file_nominal = TFile::Open("root/2016-02-26/2016-02-26_74x_template_erfcexp/scalefactors_ele_74x.root");
  TFile *dgh_file_AltBkg  = TFile::Open("root/2016-02-26/2016-02-26_74x_template_exp/scalefactors_ele_74x.root");
  TFile *dgh_file_AltSig  = TFile::Open("root/2016-02-26/2016-02-26_74x_BWCBPlusVoigt_erfcexp/scalefactors_ele_74x.root");
  TFile *dgh_file_AltMC   = TFile::Open("root/2016-02-26/2016-02-26_74x_template_erfcexp_LO/scalefactors_ele_74x.root");
  TFile *dgh_file_diffTag = TFile::Open("root/2016-02-26/2016-02-26_74x_template_erfcexp_diffTag/scalefactors_ele_74x.root");
  assert(
    dgh_file_nominal ->IsOpen() && !dgh_file_nominal ->IsZombie() &&
    dgh_file_AltBkg  ->IsOpen() && !dgh_file_AltBkg  ->IsZombie() && 
    dgh_file_AltSig  ->IsOpen() && !dgh_file_AltSig  ->IsZombie() && 
    dgh_file_AltMC   ->IsOpen() && !dgh_file_AltMC   ->IsZombie() && 
    dgh_file_diffTag ->IsOpen() && !dgh_file_diffTag ->IsZombie() 
  );
  TH2D *effData      = (TH2D*) dgh_file_nominal->Get("eff_data_Tight_ele");
  TH2D *effMC        = (TH2D*) dgh_file_nominal->Get("eff_mc_Tight_ele");
  TH2D *effAltBkg    = (TH2D*) dgh_file_AltBkg->Get("eff_data_Tight_ele");
  TH2D *effAltSig    = (TH2D*) dgh_file_AltSig->Get("eff_data_Tight_ele");
  TH2D *effAltMC     = (TH2D*) dgh_file_AltMC->Get("eff_mc_Tight_ele");
  TH2D *effAltTagSel = (TH2D*) dgh_file_diffTag->Get("eff_data_Tight_ele");
  
  TFile *egm_file = TFile::Open("root/CutBasedID_TightWP_fromTemplates_withSyst_Final.txt_SF2D.root");
  assert(egm_file->IsOpen() && !egm_file->IsZombie());
  FILE *dgh_textfile = fopen("DGH_scalefactors_ele_74x_2016-03-09.txt", "w");
  assert(dgh_textfile != NULL);
  TH2D *egm_sf = (TH2D*) egm_file->Get("EGamma_SF2D");
  for(int i = 1; i<=egm_sf->GetNbinsX(); i++) { for(int j = 1; j<=egm_sf->GetNbinsY(); j++) {
    double etacenter = egm_sf->GetXaxis()->GetBinCenter(i);
    double ptcenter = egm_sf->GetYaxis()->GetBinCenter(j);
    int dghbin = effData->FindBin( TMath::Abs(etacenter), ptcenter);
    int egmbin = egm_sf->GetBin(i,j);
    fprintf(dgh_textfile, "%-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",
      egm_sf->GetXaxis()->GetBinLowEdge(i),   //etaMin
      egm_sf->GetXaxis()->GetBinUpEdge(i),    //etaMax
      egm_sf->GetYaxis()->GetBinLowEdge(j),   //ptMin
      egm_sf->GetYaxis()->GetBinUpEdge(j),    //ptMax
      effData->GetBinContent(dghbin),         //effData
      effData->GetBinError(dghbin),           //statErrData
      effMC->GetBinContent(dghbin),           //effMC
      effMC->GetBinError(dghbin),             //statErrMC
      effAltBkg->GetBinContent(dghbin),       //effAltBkgModel
      effAltSig->GetBinContent(dghbin),       //effAltSigModel
      effAltMC->GetBinContent(dghbin),        //effAltMCSignal
      effAltTagSel->GetBinContent(dghbin),    //effAltTagSel
      -1.,                                    //effPUup
      -1.,                                    //effPUdown
      -1.                                     //effAltFitRange
    );

  }}
  assert(fclose(dgh_textfile) ==0);

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void egamma_comparison() {
  gStyle->SetOptStat(0);
  mitPalette2();
  gStyle->SetPaintTextFormat("4.3f");
  TPaletteAxis *palette_axis;
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);

  TFile *egm_file = TFile::Open("root/CutBasedID_TightWP_fromTemplates_withSyst_Final.txt_SF2D.root");
  TFile *dgh_file = TFile::Open("root/DGH_scalefactors_ele_74x_2016-03-09.root");

  TH2D *egm_sf = (TH2D*) egm_file->Get("EGamma_SF2D");
  TH2D *dgh_sf = (TH2D*) dgh_file->Get("scalefactors_Tight_ele");
  TH2D *dgh_sf_stat_hi = (TH2D*) dgh_file->Get("scalefactors_Tight_ele_stat_error_hi");
  TH2D *dgh_sf_stat_lo = (TH2D*) dgh_file->Get("scalefactors_Tight_ele_stat_error_lo");
  TH2D *dgh_sf_syst_combined = (TH2D*) dgh_file->Get("scalefactors_Tight_ele_syst_error_combined");

  TH2D *sf_differences = (TH2D*) egm_sf->Clone();
  sf_differences->Reset();
  sf_differences->SetName("sf_differences");
  sf_differences->SetTitle("Difference of tight ID/iso scale factors (DGH minus Official)");
  TH2D *error_differences = (TH2D*) egm_sf->Clone();
  error_differences->Reset();
  error_differences->SetName("error_differences");
  error_differences->SetTitle("Difference of tight ID/iso scale factor uncertainties (DGH minus Official)");

  for(int i = 1; i<=egm_sf->GetNbinsX(); i++) { for(int j = 1; j<=egm_sf->GetNbinsY(); j++) {
    double etacenter = egm_sf->GetXaxis()->GetBinCenter(i);
    double ptcenter = egm_sf->GetYaxis()->GetBinCenter(j);
    int dghbin = dgh_sf->FindBin( TMath::Abs(etacenter), ptcenter);
    int egmbin = egm_sf->GetBin(i,j);
    double dgh_error = sqrt(
      pow(TMath::Max(dgh_sf_stat_hi->GetBinContent(dghbin), dgh_sf_stat_lo->GetBinContent(dghbin)), 2) +
      pow(dgh_sf_syst_combined->GetBinContent(dghbin), 2)
    );
    double sf_difference = dgh_sf->GetBinContent(dghbin) - egm_sf->GetBinContent(egmbin);
    double error_difference = dgh_error - egm_sf->GetBinError(egmbin);
    double sf_difference_error = sqrt(
      pow(dgh_error , 2) +
      pow(egm_sf->GetBinError(egmbin), 2)
    );
    sf_differences->SetBinContent(egmbin, sf_difference);
    sf_differences->SetBinError(egmbin, sf_difference_error);
    error_differences->SetBinContent(egmbin, error_difference);

  }}
  TCanvas *c_sf_differences = new TCanvas("c_sf_differences", "Difference of scale factors", 1600, 600);
  c_sf_differences->SetLogy();
  sf_differences->GetYaxis()->SetMoreLogLabels();
  sf_differences->Draw("TEXTE COLZ");
  gPad->Update();
  sf_differences->GetXaxis()->SetTitle("#eta_{SC}");
  sf_differences->GetXaxis()->SetTitleOffset(0.9);
  sf_differences->GetXaxis()->SetTitleSize(0.04);
  sf_differences->GetXaxis()->SetLabelSize(0.02);
  sf_differences->GetYaxis()->SetTitle("p_{T} [GeV]");
  sf_differences->GetYaxis()->SetTitleOffset(0.9);
  sf_differences->GetYaxis()->SetTitleSize(0.04);
  sf_differences->GetYaxis()->SetLabelSize(0.02);
  sf_differences->GetYaxis()->SetRangeUser(10,200);
  sf_differences->SetMinimum(-.2);
  sf_differences->SetMaximum(0.2);
  sf_differences->SetMarkerSize(1.3);
  palette_axis = (TPaletteAxis*) sf_differences->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_sf_differences->Update();
  c_sf_differences->Print("sf_differences.png");

  TCanvas *c_error_differences = new TCanvas("c_error_differences", "Difference of scale factor errors", 1600, 600);
  c_error_differences->SetLogy();
  error_differences->GetYaxis()->SetMoreLogLabels();
  error_differences->Draw("TEXT COLZ");
  gPad->Update();
  error_differences->GetXaxis()->SetTitle("#eta_{SC}");
  error_differences->GetXaxis()->SetTitleOffset(0.9);
  error_differences->GetXaxis()->SetTitleSize(0.04);
  error_differences->GetXaxis()->SetLabelSize(0.02);
  error_differences->GetYaxis()->SetTitle("p_{T} [GeV]");
  error_differences->GetYaxis()->SetTitleOffset(0.9);
  error_differences->GetYaxis()->SetTitleSize(0.04);
  error_differences->GetYaxis()->SetLabelSize(0.02);
  error_differences->GetYaxis()->SetRangeUser(10,200);
  error_differences->SetMinimum(-.2);
  error_differences->SetMaximum(0.2);
  error_differences->SetMarkerSize(1.3);
  palette_axis = (TPaletteAxis*) error_differences->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_error_differences->Update();
  c_error_differences->Print("error_differences.png");

}
