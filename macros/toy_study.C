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
void compare_methods(string toys_dir, string plots_dir, string method1, string method2, const unsigned int numtoys, string toy_dir_format="toy%06d") {
  // This function looks for folders like toy012345 in <toys_dir>
  // then compares the efficiencies in files with name formatted like eff_<method1>_etapt_<bin>.txt
  // Those file contents are as such: <efficiency> <error-high> <error-low>
  // The output is plot image files in plots_dir.
  //
  // Currently the eta and pT binnings are hardcoded.
  // Other options besides eta-pT phase space bins will be added later.

  gStyle->SetOptStat(2210);
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);
  if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
  if( toys_dir[toys_dir.size()-1]  != '/' )   toys_dir  = toys_dir + "/";
  TH2D *syst_from_Ztest = new TH2D("syst_from_Ztest", "#LT (#varepsilon_{2} #minus #varepsilon_{1})/#sqrt{#sigma_{1}^{2}+#sigma_{2}^{2}} #GT_{toys} #LT #sqrt{#sigma_{1}^{2}+#sigma_{2}^{2}} #GT_{toys} #cbar_{(|#eta|, p_{T})}", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);

  char title[512], filename[1024];
  for(int j=0; j<n_ele_pt_bins; j++) { for(int i=0; i<n_ele_eta_bins; i++) { 
    int ibin = i + n_ele_eta_bins*j;
    TGraphAsymmErrors *toy_eff_method1=new TGraphAsymmErrors(numtoys);
    TGraphAsymmErrors *toy_eff_method2=new TGraphAsymmErrors(numtoys);
    sprintf(title, "dN/d#varepsilon_{1} #cbar_{(|#eta|, p_{T}) = (%.4f, %d)}^{(%.4f, %d)}", ele_eta_bins[i], (int) ele_pt_bins[j], ele_eta_bins[i+1], (int) ele_pt_bins[j+1]);
    TH1D *dNdE_method1 = new TH1D("dNdE_method1", title, 1000,0,1);
    sprintf(title, "dN/d#varepsilon_{2} #cbar_{(|#eta|, p_{T}) = (%.4f, %d)}^{(%.4f, %d)}", ele_eta_bins[i], (int) ele_pt_bins[j], ele_eta_bins[i+1], (int) ele_pt_bins[j+1]);
    TH1D *dNdE_method2 = new TH1D("dNdE_method2", title, 1000,0,1);
    sprintf(title, "Z (%.4f < |#eta| < %.4f, %d < p_{T} < %d)", ele_eta_bins[i], ele_eta_bins[i+1], (int) ele_pt_bins[j], (int) ele_pt_bins[j+1]);
    TH1D *eff_Ztest = new TH1D("eff_Ztest", title, 400, -20., 20.); // +/- 20 sigma
    sprintf(title, "Distribution of uncertainties for Z-test #sqrt{#sigma_{1}^{2}+#sigma_{2}^{2}} (%d toys)", numtoys);
    TH1D *eff_sigma = new TH1D("eff_sigma", title, 50, 0, 0.1);
    for(unsigned int toynum=0; toynum<numtoys; toynum++) {
      double eff1, errh1, errl1, eff2, errh2, errl2, sigma=.0001, Ztest=0;

      sprintf(filename, (toys_dir+toy_dir_format+"/eff_%s_etapt_%d.txt").c_str(), toynum, method1.c_str(), ibin);
      //printf("trying to open \"%s\"\n", filename);
      std::ifstream file_method1(filename);
      if(!file_method1.is_open()) continue;
      file_method1 >> eff1 >> errh1 >> errl1;
      toy_eff_method1->SetPoint(toynum, toynum, eff1);
      toy_eff_method1->SetPointError(toynum, 0,0, errl1, errh1);

      sprintf(filename, (toys_dir+toy_dir_format+"/eff_%s_etapt_%d.txt").c_str(), toynum, method2.c_str(), ibin);
      //printf("trying to open \"%s\"\n", filename);
      std::ifstream file_method2(filename);
      if(!file_method2.is_open()) continue;
      file_method2 >> eff2 >> errh2 >> errl2;
      toy_eff_method2->SetPoint(toynum, toynum, eff2);
      toy_eff_method2->SetPointError(toynum, 0,0, errl2, errh2);

      if(eff2<eff1)      sigma=sqrt(errl1*errl1 + errh2*errh2);
      else if(eff2>eff1) sigma=sqrt(errh1*errh1 + errl2*errl2);
      else printf("wow amazing!!\n");
      Ztest=(eff2-eff1)/sigma;

      dNdE_method1->Fill(eff1);
      dNdE_method2->Fill(eff2);
      eff_Ztest->Fill(Ztest);
      eff_sigma->Fill(sigma);
    }
    syst_from_Ztest->SetBinContent(syst_from_Ztest->GetBin(i+1,j+1), TMath::Abs(eff_Ztest->GetMean() * eff_sigma->GetMean()));
    syst_from_Ztest->SetBinError(syst_from_Ztest->GetBin(i+1,j+1), TMath::Abs(eff_Ztest->GetMeanError() * eff_sigma->GetMeanError()));

    // start making plots
    gStyle->SetOptStat(2210);
    TCanvas *canvas = new TCanvas;
    canvas->SetTopMargin(.15);
    dNdE_method1->Draw();
    dNdE_method1->GetXaxis()->SetRangeUser(
      dNdE_method1->GetMean() - 5.*dNdE_method1->GetRMS(),
      dNdE_method1->GetMean() + 5.*dNdE_method1->GetRMS()
    );
    dNdE_method1->GetXaxis()->SetTitle("d#varepsilon_{1}");
    dNdE_method1->GetYaxis()->SetTitle("N / 0.001");
    canvas->Update();
    sprintf(filename, (plots_dir+"toystudy_dNdE_%s_etapt_%d.png").c_str(), method1.c_str(), ibin);
    canvas->Print(filename);
    
    dNdE_method2->Draw();
    dNdE_method2->GetXaxis()->SetRangeUser(
      dNdE_method2->GetMean() - 5.*dNdE_method2->GetRMS(),
      dNdE_method2->GetMean() + 5.*dNdE_method2->GetRMS()
    );
    dNdE_method2->GetXaxis()->SetTitle("d#varepsilon_{2}");
    dNdE_method2->GetYaxis()->SetTitle("N / 0.001");
    canvas->Update();
    sprintf(filename, (plots_dir+"toystudy_dNdE_%s_etapt_%d.png").c_str(), method2.c_str(), ibin);
    canvas->Print(filename);
    
    sprintf(filename, (plots_dir+"toystudy_eff_Ztest_%s_versus_%s_etapt_%d.png").c_str(), method1.c_str(), method2.c_str(), ibin);
    eff_Ztest->Draw();    canvas->Print(filename);
    
    sprintf(filename, (plots_dir+"toystudy_eff_sigma_%s_versus_%s_etapt_%d.png").c_str(), method1.c_str(), method2.c_str(), ibin);
    eff_sigma->Draw();    canvas->Print(filename);

    TCanvas *c_toy_eff = new TCanvas("c_toy_eff", "c_toy_eff", 1600,600);
    c_toy_eff->SetRightMargin(.05);
    c_toy_eff->SetLeftMargin(.05);
    sprintf(filename, (plots_dir+"toystudy_efficiency_%s_versus_%s_etapt_%d.png").c_str(), method1.c_str(), method2.c_str(), ibin);
    toy_eff_method2->SetLineColor(mit_red);
    toy_eff_method2->SetFillColorAlpha(mit_red,.4);
    toy_eff_method2->SetFillStyle(1001);

    toy_eff_method1->SetLineColor(1);
    toy_eff_method1->SetFillColorAlpha(mit_gray,.4);
    toy_eff_method1->SetFillStyle(1001);

    sprintf(title, "Toy efficiency comparison (%.4f < |#eta| < %.4f, %d < p_{T} < %d)", ele_eta_bins[i], ele_eta_bins[i+1], (int) ele_pt_bins[j], (int) ele_pt_bins[j+1]);
    toy_eff_method2->SetTitle(title);
    toy_eff_method2->Draw("AL E3");
    toy_eff_method2->GetXaxis()->SetRangeUser(0, numtoys);
    toy_eff_method2->GetXaxis()->SetTitle("n_{toy}");
    toy_eff_method2->GetYaxis()->SetRangeUser(
      TMath::Min( dNdE_method1->GetMean() - 4.*dNdE_method1->GetRMS(), dNdE_method2->GetMean() - 4.*dNdE_method2->GetRMS()), 
      TMath::Max( dNdE_method1->GetMean() + 6.*dNdE_method1->GetRMS(), dNdE_method2->GetMean() + 6.*dNdE_method2->GetRMS())
    );
    toy_eff_method1->Draw("L SAME E3");
    TLegend *l_eff = new TLegend(.7,.7,.9,.85);
    l_eff->AddEntry(toy_eff_method1, method1.c_str(), "lf");
    l_eff->AddEntry(toy_eff_method2, method2.c_str(), "lf");
    l_eff->SetFillStyle(1001);
    l_eff->SetFillColorAlpha(0, .3);
    l_eff->Draw("SAME");
    c_toy_eff->Print(filename);

    delete dNdE_method1;
    delete dNdE_method2;
    delete eff_Ztest;
    delete eff_sigma;
    delete toy_eff_method1;
    delete toy_eff_method2;
    delete canvas;
    delete c_toy_eff;

  }}

  TCanvas *c_syst_from_Ztest = new TCanvas("c_syst_from_Ztest", "c_syst_from_Ztest", 800, 800);
  c_syst_from_Ztest->SetRightMargin(.15);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  mitPalette2();
  TPaletteAxis *palette_axis;
  syst_from_Ztest->Draw("TEXT COLZ");
  gPad->Update();
  syst_from_Ztest->GetXaxis()->SetTitle("| #eta |");
  syst_from_Ztest->GetXaxis()->SetTitleOffset(0.9);
  syst_from_Ztest->GetXaxis()->SetTitleSize(0.04);
  syst_from_Ztest->GetXaxis()->SetLabelSize(0.02);
  syst_from_Ztest->GetYaxis()->SetTitle("p_{T} [GeV]");
  syst_from_Ztest->GetYaxis()->SetTitleOffset(0.9);
  syst_from_Ztest->GetYaxis()->SetTitleSize(0.04);
  syst_from_Ztest->GetYaxis()->SetLabelSize(0.02);
  syst_from_Ztest->GetYaxis()->SetRangeUser(10,100);
  syst_from_Ztest->SetMinimum(-.2);
  syst_from_Ztest->SetMaximum(0.2);
  syst_from_Ztest->SetMarkerSize(1.3);
  palette_axis = (TPaletteAxis*) syst_from_Ztest->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_syst_from_Ztest->Update();
  sprintf(filename, (plots_dir+"toystudy_syst_from_Ztest_%s_versus_%s.png").c_str(), method1.c_str(), method2.c_str());
  c_syst_from_Ztest->Print(filename);
  delete syst_from_Ztest;
  delete c_syst_from_Ztest;

}
void check_toy() {
  printf("wow\n");
}
