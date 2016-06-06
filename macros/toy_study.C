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
Float_t mu_pt_bins[] = {10,20,25,30,40,50,60,120};
Float_t mu_eta_bins[] = {0, 0.9, 1.2, 2.1, 2.4};
Int_t n_mu_pt_bins=7;
Int_t n_mu_eta_bins=4;
int marker_colors[] = {1, mit_gray, mit_red, 97, 91, 8, 60};
int marker_styles[] = {20, 21, 22, 23, 33, 34, 20};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compare_methods(string toys_dir, string plots_dir, string root_dir, string method1, string method2, const unsigned int numtoys, string toy_dir_format="toy%06d", string flavor="electrons") {
  // This function looks for folders like toy012345 in <toys_dir>
  // then compares the efficiencies in files with name formatted like eff_<method1>_etapt_<bin>.txt
  // Those file contents are as such: <efficiency> <error-high> <error-low>
  // The output is plot image files in plots_dir.
  //
  // Currently the eta and pT binnings are hardcoded.
  // Other options besides eta-pT phase space bins will be added later.

  if(flavor != "electrons" && flavor != "muons") return;
  TColor *col_mit_red  = new TColor(mit_red,  163/255., 31/255.,  52/255.);
  TColor *col_mit_gray = new TColor(mit_gray, 138/255., 139/255., 140/255.);
  if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
  if( root_dir[root_dir.size()-1]  != '/' )   root_dir  = root_dir + "/";
  if( toys_dir[toys_dir.size()-1]  != '/' )   toys_dir  = toys_dir + "/";
  system(("mkdir -p "+plots_dir).c_str());
  system(("mkdir -p "+root_dir).c_str());
  TH2D *syst_from_Ztest, *mean_eff_method1, *mean_eff_method2;
  if(flavor=="electrons") {
    syst_from_Ztest = new TH2D("syst_from_Ztest", "#LT (#varepsilon_{2} #minus #varepsilon_{1})/#sqrt{#sigma_{1}^{2}+#sigma_{2}^{2}} #GT_{toys} #LT #sqrt{#sigma_{1}^{2}+#sigma_{2}^{2}} #GT_{toys} #cbar_{(|#eta|, p_{T})}", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
    mean_eff_method1 = new TH2D("mean_eff_method1", "#LT #varepsilon_{1} #GT_{toys}", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
    mean_eff_method2 = new TH2D("mean_eff_method2", "#LT #varepsilon_{2} #GT_{toys}", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  } else if(flavor=="muons") {
    syst_from_Ztest = new TH2D("syst_from_Ztest", "#LT (#varepsilon_{2} #minus #varepsilon_{1})/#sqrt{#sigma_{1}^{2}+#sigma_{2}^{2}} #GT_{toys} #LT #sqrt{#sigma_{1}^{2}+#sigma_{2}^{2}} #GT_{toys} #cbar_{(|#eta|, p_{T})}", n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
    mean_eff_method1 = new TH2D("mean_eff_method1", "#LT #varepsilon_{1} #GT_{toys}", n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
    mean_eff_method2 = new TH2D("mean_eff_method2", "#LT #varepsilon_{2} #GT_{toys}", n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  }

  char title[512], filename[1024];
  bool missing_files_method1[numtoys], missing_files_method2[numtoys];
  for(unsigned int toynum=0; toynum<numtoys; toynum++) {
    missing_files_method1[toynum]=false;
    missing_files_method2[toynum]=false;
  }
  for(int j=1; j<=mean_eff_method1->GetNbinsY(); j++) { for(int i=1; i<=mean_eff_method1->GetNbinsX(); i++) { 
    int ibin = (i-1) + mean_eff_method1->GetNbinsX()*(j-1);
    double eta_min = mean_eff_method1->GetXaxis()->GetBinLowEdge(i);
    double eta_max = mean_eff_method1->GetXaxis()->GetBinUpEdge(i);
    double pt_min = mean_eff_method1->GetYaxis()->GetBinLowEdge(j);
    double pt_max = mean_eff_method1->GetYaxis()->GetBinUpEdge(j);
    TGraphAsymmErrors *toy_eff_method1=new TGraphAsymmErrors(numtoys);
    TGraphAsymmErrors *toy_eff_method2=new TGraphAsymmErrors(numtoys);
    sprintf(title, "dN/d#varepsilon_{1} #cbar_{(|#eta|, p_{T}) = (%.4f, %d)}^{(%.4f, %d)}", eta_min, (int) pt_min, eta_max, (int) pt_max);
    TH1D *dNdE_method1 = new TH1D("dNdE_method1", title, 1000,0,1);
    sprintf(title, "dN/d#varepsilon_{2} #cbar_{(|#eta|, p_{T}) = (%.4f, %d)}^{(%.4f, %d)}", eta_min, (int) pt_min, eta_max, (int) pt_max);
    TH1D *dNdE_method2 = new TH1D("dNdE_method2", title, 1000,0,1);
    sprintf(title, "Z (%.4f < |#eta| < %.4f, %d < p_{T} < %d)", eta_min, eta_max, (int) pt_min, (int) pt_max);
    TH1D *eff_Ztest = new TH1D("eff_Ztest", title, 400, -20., 20.); // +/- 20 sigma
    sprintf(title, "Distribution of uncertainties for Z-test #sqrt{#sigma_{1}^{2}+#sigma_{2}^{2}} (%d toys)", numtoys);
    TH1D *eff_sigma = new TH1D("eff_sigma", title, 50, 0, 0.1);
    for(unsigned int toynum=0; toynum<numtoys; toynum++) {
      double eff1, errh1, errl1, eff2, errh2, errl2, sigma=.0001, Ztest=0;

      sprintf(filename, (toys_dir+toy_dir_format+"/eff_%s_etapt_%d.txt").c_str(), toynum, method1.c_str(), ibin);
      //printf("trying to open \"%s\"\n", filename);
      std::ifstream file_method1(filename);
      if(!file_method1.is_open()) { 
        missing_files_method1[toynum]=true;
        printf("ERROR: File is missing! (method %s, toy #%d, bin #%d)\n", method1.c_str(), toynum, ibin);
        continue;
      }
      //assert(file_method1.is_open());
      file_method1 >> eff1 >> errl1 >> errh1;
      toy_eff_method1->SetPoint(toynum, toynum, eff1);
      toy_eff_method1->SetPointError(toynum, 0,0, errl1, errh1);

      sprintf(filename, (toys_dir+toy_dir_format+"/eff_%s_etapt_%d.txt").c_str(), toynum, method2.c_str(), ibin);
      //printf("trying to open \"%s\"\n", filename);
      std::ifstream file_method2(filename);
      if(!file_method2.is_open()) { 
        missing_files_method2[toynum]=true;
        printf("ERROR: File is missing! (method %s, toy #%d, bin #%d)\n", method2.c_str(), toynum, ibin);
        continue;
      }
      //assert(file_method2.is_open());
      file_method2 >> eff2 >> errl2 >> errh2;
      toy_eff_method2->SetPoint(toynum, toynum, eff2);
      toy_eff_method2->SetPointError(toynum, 0,0, errl2, errh2);

      if(eff2<eff1)      sigma=sqrt(errl1*errl1 + errh2*errh2);
      else if(eff2>eff1) sigma=sqrt(errh1*errh1 + errl2*errl2);
      //else printf("wow amazing!!\n");
      Ztest=(eff2-eff1)/sigma;

      dNdE_method1->Fill(eff1);
      dNdE_method2->Fill(eff2);
      eff_Ztest->Fill(Ztest);
      eff_sigma->Fill(sigma);
    }
    syst_from_Ztest->SetBinContent(syst_from_Ztest->GetBin(i,j), TMath::Abs(eff_Ztest->GetMean() * eff_sigma->GetMean()));
    syst_from_Ztest->SetBinError(syst_from_Ztest->GetBin(i,j), TMath::Abs(eff_Ztest->GetMeanError() * eff_sigma->GetMeanError()));
    mean_eff_method1->SetBinContent( mean_eff_method1->GetBin(i,j), dNdE_method1->GetMean() );
    mean_eff_method1->SetBinError( mean_eff_method1->GetBin(i,j), dNdE_method1->GetRMS() );
    mean_eff_method2->SetBinContent( mean_eff_method2->GetBin(i,j), dNdE_method2->GetMean() );
    mean_eff_method2->SetBinError( mean_eff_method2->GetBin(i,j), dNdE_method2->GetRMS() );

    // start making plots
    gStyle->SetOptStat(2210);
    TCanvas *canvas = new TCanvas;
    canvas->UseCurrentStyle();
    canvas->SetTopMargin(.15);
    dNdE_method1->SetFillStyle(3004);
    dNdE_method1->SetLineColor(1);
    dNdE_method1->SetFillColor(1);
    dNdE_method2->SetFillStyle(3005);
    dNdE_method2->SetLineColor(mit_red);
    dNdE_method2->SetFillColor(mit_red);
    dNdE_method1->Draw();
    dNdE_method1->GetXaxis()->SetRangeUser(
      TMath::Min(dNdE_method1->GetMean() - 5.*dNdE_method1->GetRMS(), dNdE_method2->GetMean() - 5.*dNdE_method2->GetRMS()),
      TMath::Max(dNdE_method1->GetMean() + 5.*dNdE_method1->GetRMS(), dNdE_method2->GetMean() + 5.*dNdE_method2->GetRMS())
    );
    dNdE_method1->GetXaxis()->SetTitle("d#varepsilon");
    sprintf(title, "N / %f", 1./numtoys);
    dNdE_method1->GetYaxis()->SetTitle(title);
    canvas->Update();
    sprintf(filename, (plots_dir+"toystudy_dNdE_%s_etapt_%d.png").c_str(), method1.c_str(), ibin);
    canvas->Print(filename);
    
    dNdE_method2->Draw();
    dNdE_method2->GetXaxis()->SetTitle("d#varepsilon");
    sprintf(title, "N / %f", 1./numtoys);
    dNdE_method2->GetYaxis()->SetTitle(title);
    dNdE_method2->GetXaxis()->SetRangeUser(
      TMath::Min(dNdE_method2->GetMean() - 5.*dNdE_method2->GetRMS(), dNdE_method2->GetMean() - 5.*dNdE_method2->GetRMS()),
      TMath::Max(dNdE_method2->GetMean() + 5.*dNdE_method2->GetRMS(), dNdE_method2->GetMean() + 5.*dNdE_method2->GetRMS())
    );
    canvas->Update();
    sprintf(filename, (plots_dir+"toystudy_dNdE_%s_etapt_%d.png").c_str(), method2.c_str(), ibin);
    canvas->Print(filename);
    
    sprintf(title, "dN/d#varepsilon comparison (%.4f < |#eta| < %.4f, %d < p_{T} < %d)", eta_min, eta_max, (int) pt_min, (int) pt_max);
    dNdE_method1->SetTitle(title);
    
    gStyle->SetOptStat(0);
    canvas->UseCurrentStyle();
    dNdE_method1->Draw();
    dNdE_method2->Draw("SAME");
    TLegend *l_dNdE = new TLegend(.65,.65,.85,.8);
    l_dNdE->AddEntry(dNdE_method1, method1.c_str(), "lf");
    l_dNdE->AddEntry(dNdE_method2, method2.c_str(), "lf");
    l_dNdE->SetFillStyle(1001);
    l_dNdE->SetFillColorAlpha(0, .3);
    l_dNdE->Draw("SAME");
    sprintf(filename, (plots_dir+"toystudy_dNdE_%s_versus_%s_etapt_%d.png").c_str(), method1.c_str(), method2.c_str(), ibin);
    canvas->Print(filename);
    
    gStyle->SetOptStat(2210);
    canvas->UseCurrentStyle();
    sprintf(filename, (plots_dir+"toystudy_eff_Ztest_%s_versus_%s_etapt_%d.png").c_str(), method1.c_str(), method2.c_str(), ibin);
    eff_Ztest->Draw();    canvas->Print(filename);
    
    sprintf(filename, (plots_dir+"toystudy_eff_sigma_%s_versus_%s_etapt_%d.png").c_str(), method1.c_str(), method2.c_str(), ibin);
    eff_sigma->Draw();    canvas->Print(filename);

    gStyle->SetOptStat(0);
    TCanvas *c_toy_eff = new TCanvas("c_toy_eff", "c_toy_eff", 1600,600);
    canvas->UseCurrentStyle();
    c_toy_eff->SetRightMargin(.05);
    c_toy_eff->SetLeftMargin(.05);
    sprintf(filename, (plots_dir+"toystudy_efficiency_%s_versus_%s_etapt_%d.png").c_str(), method1.c_str(), method2.c_str(), ibin);
    toy_eff_method2->SetLineColor(mit_red);
    toy_eff_method2->SetFillColorAlpha(mit_red,.4);
    toy_eff_method2->SetFillStyle(1001);

    toy_eff_method1->SetLineColor(1);
    toy_eff_method1->SetFillColorAlpha(mit_gray,.4);
    toy_eff_method1->SetFillStyle(1001);

    sprintf(title, "Toy efficiency comparison (%.4f < |#eta| < %.4f, %d < p_{T} < %d)", eta_min, eta_max, (int) pt_min, (int) pt_max);
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
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  TPaletteAxis *palette_axis;

  mitPalette();
  TCanvas *c_mean_eff_method1 = new TCanvas("c_mean_eff_method1", "c_mean_eff_method1", 800, 800);
  c_mean_eff_method1->SetRightMargin(.15);
  mean_eff_method1->Draw("TEXTE COLZ");
  gPad->Update();
  mean_eff_method1->GetXaxis()->SetTitle("| #eta |");
  mean_eff_method1->GetXaxis()->SetTitleOffset(0.9);
  mean_eff_method1->GetXaxis()->SetTitleSize(0.04);
  mean_eff_method1->GetXaxis()->SetLabelSize(0.02);
  mean_eff_method1->GetYaxis()->SetTitle("p_{T} [GeV]");
  mean_eff_method1->GetYaxis()->SetTitleOffset(0.9);
  mean_eff_method1->GetYaxis()->SetTitleSize(0.04);
  mean_eff_method1->GetYaxis()->SetLabelSize(0.02);
  mean_eff_method1->GetYaxis()->SetRangeUser(10,100);
  mean_eff_method1->SetMinimum(0);
  mean_eff_method1->SetMaximum(1);
  mean_eff_method1->SetMarkerSize(1.3);
  palette_axis = (TPaletteAxis*) mean_eff_method1->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_mean_eff_method1->Update();
  sprintf(filename, (plots_dir+"toystudy_mean_eff_%s.png").c_str(), method1.c_str());
  c_mean_eff_method1->Print(filename);
  
  TCanvas *c_mean_eff_method2 = new TCanvas("c_mean_eff_method2", "c_mean_eff_method2", 800, 800);
  c_mean_eff_method2->SetRightMargin(.15);
  mean_eff_method2->Draw("TEXTE COLZ");
  gPad->Update();
  mean_eff_method2->GetXaxis()->SetTitle("| #eta |");
  mean_eff_method2->GetXaxis()->SetTitleOffset(0.9);
  mean_eff_method2->GetXaxis()->SetTitleSize(0.04);
  mean_eff_method2->GetXaxis()->SetLabelSize(0.02);
  mean_eff_method2->GetYaxis()->SetTitle("p_{T} [GeV]");
  mean_eff_method2->GetYaxis()->SetTitleOffset(0.9);
  mean_eff_method2->GetYaxis()->SetTitleSize(0.04);
  mean_eff_method2->GetYaxis()->SetLabelSize(0.02);
  mean_eff_method2->GetYaxis()->SetRangeUser(10,100);
  mean_eff_method2->SetMinimum(0);
  mean_eff_method2->SetMaximum(1);
  mean_eff_method2->SetMarkerSize(1.3);
  palette_axis = (TPaletteAxis*) mean_eff_method2->GetListOfFunctions()->FindObject("palette"); 
  palette_axis->SetLabelSize(0.02);
  c_mean_eff_method2->Update();
  sprintf(filename, (plots_dir+"toystudy_mean_eff_%s.png").c_str(), method2.c_str());
  c_mean_eff_method2->Print(filename);

  TCanvas *c_syst_from_Ztest = new TCanvas("c_syst_from_Ztest", "c_syst_from_Ztest", 800, 800);
  c_syst_from_Ztest->SetRightMargin(.15);
  mitPalette2();
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
  sprintf(filename, (root_dir+"toystudy_%s_versus_%s.root").c_str(), method1.c_str(), method2.c_str());

  TFile *file = TFile::Open(filename, "RECREATE");
  syst_from_Ztest->Write();
  mean_eff_method1->Write();
  mean_eff_method2->Write();
  file->Close();

  delete syst_from_Ztest;
  delete c_syst_from_Ztest;
  delete mean_eff_method1;
  delete c_mean_eff_method1;
  delete mean_eff_method2;
  delete c_mean_eff_method2;

  vector<unsigned int> toys_to_redo_method1, toys_to_redo_method2;
  for(unsigned int toynum=0; toynum<numtoys; toynum++) {
    if(missing_files_method1[toynum]) toys_to_redo_method1.push_back(toynum);
    if(missing_files_method2[toynum]) toys_to_redo_method2.push_back(toynum);
  }
  if(toys_to_redo_method1.size()>0 || toys_to_redo_method2.size()>0) {
    printf("\nSUMMARY: Some files were missing!\n");
    if(toys_to_redo_method1.size()>0) { 
      printf("Files missing for method %s:\n", method1.c_str());
      for(unsigned int i=0; i<(unsigned int)toys_to_redo_method1.size(); i++) {
        printf("%d ", toys_to_redo_method1[i]);
      }
      printf("\n");
    }
    if(toys_to_redo_method2.size()>0) { 
      printf("Files missing for method %s:\n", method2.c_str());
      for(unsigned int i=0; i<(unsigned int)toys_to_redo_method2.size(); i++) {
        printf("%d ", toys_to_redo_method2[i]);
      }
      printf("\n");
    }
  }

}
void combined_toy_systematics(string sf_filename, string plots_dir, string root_dir, string flavor="electrons") {
  assert(flavor=="electrons" || flavor=="muons");
  unsigned int n_working_points;
  char *flavors_short, *flavors_caps, *flavors_lowercase, **selections_caps, **selections_lowercase;
  if(flavor=="electrons") {
    n_working_points=4;
    flavors_short = "ele";
    flavors_caps = "Electron";
    flavors_lowercase = "electron";
    selections_caps      = new char*[n_working_points];
    selections_lowercase = new char*[n_working_points];
    selections_caps     [0] = "Veto"  ;
    selections_caps     [1] = "Loose" ;
    selections_caps     [2] = "Medium";
    selections_caps     [3] = "Tight" ;
    selections_lowercase[0] = "veto"  ;
    selections_lowercase[1] = "loose" ;
    selections_lowercase[2] = "medium";
    selections_lowercase[3] = "tight" ;
  } else if(flavor=="muons") {
    n_working_points=3;
    flavors_short = "mu";
    flavors_caps = "Muon";
    flavors_lowercase = "muon";
    selections_caps      = new char*[n_working_points];
    selections_lowercase = new char*[n_working_points];
    selections_caps     [0] = "Loose" ;
    selections_caps     [1] = "Medium";
    selections_caps     [2] = "Tight" ;
    selections_lowercase[0] = "loose" ;
    selections_lowercase[1] = "medium";
    selections_lowercase[2] = "tight" ;
  }
  //char const *flavors_short[]        = {"ele"       , "ele"       , "ele"       , "ele"      };
  //char const *flavors_caps[]         = {"Electron"  , "Electron"  , "Electron"  , "Electron" };
  //char const *flavors_lowercase[]    = {"electron"  , "electron"  , "electron"  , "electron" };
  //char const *selections_caps[]      = {"Veto"      , "Loose"     , "Medium"    , "Tight"    };
  //char const *selections_lowercase[] = {"veto"      , "loose"     , "medium"    , "tight"    };
  char const directory_format[]="Single%s_BaselineTo%s%s_%sTnP";
  char const diff_tag_str[]="_diffTag";

  if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
  if( root_dir[root_dir.size()-1]  != '/' ) root_dir = root_dir + "/";
  TFile *sf_file = TFile::Open(sf_filename.c_str(), "UPDATE");
  TFile *tag_cuts_file, *background_file, *signal_file, *generator_file;
  TH2D *eff_mc, *eff_mc_LO, *eff_mc_alt_tag, *eff_data;
  TFile *file_eff_mc_LO, *file_eff_mc_alt_tag;
  char filename[512], object_name[128], title[512];
  for(unsigned int iwp=0; iwp<n_working_points; iwp++) {
    if(flavor=="electrons") {
      sprintf(filename, "~/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_LO_BaselineTo%s_%sTnP/eff.root", selections_caps[iwp], flavors_lowercase);
      file_eff_mc_LO = TFile::Open(filename);
      eff_mc_LO = (TH2D*) file_eff_mc_LO->Get("hEffEtaPt");
      sprintf(filename, "~/TagAndProbe/2016-05-02_76x_electrons/template_erfcexp/DYJetsToLL_BaselineTo%s_diffTag_%sTnP/eff.root", selections_caps[iwp], flavors_lowercase);
      file_eff_mc_alt_tag=TFile::Open(filename);
      eff_mc_alt_tag = (TH2D*)file_eff_mc_alt_tag->Get("hEffEtaPt");
    } else if(flavor=="muons") {
      sprintf(filename, "~/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/DYJetsToLL_LO_BaselineTo%s_%sTnP/eff.root", selections_caps[iwp], flavors_lowercase);
      file_eff_mc_LO = TFile::Open(filename);
      eff_mc_LO = (TH2D*) file_eff_mc_LO->Get("hEffEtaPt");
      sprintf(filename, "~/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/DYJetsToLL_BaselineTo%s_diffTag_%sTnP/eff.root", selections_caps[iwp], flavors_lowercase);
      file_eff_mc_alt_tag=TFile::Open(filename);
      eff_mc_alt_tag = (TH2D*)file_eff_mc_alt_tag->Get("hEffEtaPt");
    }
    sprintf(object_name, "eff_mc_%s_%s", selections_caps[iwp], flavors_caps); 
    eff_mc = (TH2D*) sf_file->Get(object_name);
    sprintf(object_name, "eff_data_%s_%s", selections_caps[iwp], flavors_caps);
    eff_data = (TH2D*) sf_file->Get(object_name);
  
    sprintf(filename, "%sSingle%s_BaselineTo%s_diffTag_%sTnP/toystudy_template_erfcexp_diffTag_versus_template_erfcexp_diffTag.root", root_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    tag_cuts_file   = TFile::Open(filename);
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/toystudy_template_erfcexp_versus_template_exp.root", root_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    background_file = TFile::Open(filename);
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/toystudy_template_erfcexp_versus_BWCBPlusVoigt_erfcexp.root", root_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    signal_file     = TFile::Open(filename);
    //sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/toystudy_template_erfcexp_versus_template_erfcexp_LO.root", root_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    //generator_file  = TFile::Open(filename);
  
    TH2D *eff_syst_background_shape = (TH2D*) ((TH2D*) background_file->Get("syst_from_Ztest"))->Clone(object_name);
    TH2D *eff_syst_signal_shape     = (TH2D*) ((TH2D*) signal_file->Get("syst_from_Ztest"))->Clone(object_name);
    TH2D *eff_syst_generator        = (TH2D*) ((TH2D*) signal_file->Get("syst_from_Ztest"))->Clone(object_name);

    TH2D *eff_normal_tag = (TH2D*) signal_file->Get("mean_eff_method1");
    TH2D *eff_alt_tag    = (TH2D*) tag_cuts_file->Get("mean_eff_method1");

    // Create new TH2's for the SF uncertainties since Root doesn't properly handle TPaletteAxis object for cloned TH2
    unsigned int n_eta_bins = (unsigned int) eff_normal_tag->GetNbinsX();
    unsigned int n_pt_bins  = (unsigned int) eff_normal_tag->GetNbinsY();
    double eta_min = eff_normal_tag->GetXaxis()->GetBinLowEdge(1);
    double eta_max = eff_normal_tag->GetXaxis()->GetBinUpEdge(n_eta_bins);
    double pt_min  = eff_normal_tag->GetYaxis()->GetBinLowEdge(1);
    double pt_max  = eff_normal_tag->GetYaxis()->GetBinUpEdge(n_pt_bins);

    TH2D *eff_syst_tag_cuts  = (TH2D*) eff_normal_tag->Clone(object_name);
    eff_syst_tag_cuts->Add(eff_alt_tag, -1);
    
    sprintf(object_name, "scalefactors_%s_%s_syst_error_combined", selections_caps[iwp], flavors_caps);
    TH2D *sf_syst_combined = (TH2D*) eff_mc_LO->Clone(object_name);
    sprintf(object_name, "scalefactors_%s_%s_syst_error_background_shape", selections_caps[iwp], flavors_caps);
    TH2D *sf_syst_background_shape = (TH2D*) eff_mc_LO->Clone(object_name); 
    sprintf(object_name, "scalefactors_%s_%s_syst_error_signal_shape", selections_caps[iwp], flavors_caps);
    TH2D *sf_syst_signal_shape = (TH2D*) eff_mc_LO->Clone(object_name); 
    sprintf(object_name, "scalefactors_%s_%s_syst_error_generator", selections_caps[iwp], flavors_caps);
    TH2D *sf_syst_generator = (TH2D*) eff_mc_LO->Clone(object_name); 
    sprintf(object_name, "scalefactors_%s_%s_syst_error_tag_cuts", selections_caps[iwp], flavors_caps);
    TH2D *sf_syst_tag_cuts = (TH2D*) eff_mc_LO->Clone(object_name); 
  
    for(unsigned int i = 1; i <= (unsigned int) eff_syst_tag_cuts->GetNbinsX(); i++) { for(unsigned int j = 1; j <= (unsigned int) eff_syst_tag_cuts->GetNbinsY(); j++) {
      int nbin = eff_syst_tag_cuts->GetBin(i,j);
      double eff_data_value    = eff_data->GetBinContent(eff_data->FindBin(eff_syst_tag_cuts->GetXaxis()->GetBinCenter(i), eff_syst_tag_cuts->GetYaxis()->GetBinCenter(j)));
      double eff_mc_value    = eff_mc->GetBinContent(eff_mc->FindBin(eff_syst_tag_cuts->GetXaxis()->GetBinCenter(i), eff_syst_tag_cuts->GetYaxis()->GetBinCenter(j)));
      double eff_mc_alt_tag_value    = eff_mc_alt_tag->GetBinContent(eff_mc_alt_tag->FindBin(eff_syst_tag_cuts->GetXaxis()->GetBinCenter(i), eff_syst_tag_cuts->GetYaxis()->GetBinCenter(j)));
      double eff_mc_value_LO = eff_mc_LO->GetBinContent(eff_mc_LO->FindBin(eff_syst_tag_cuts->GetXaxis()->GetBinCenter(i), eff_syst_tag_cuts->GetYaxis()->GetBinCenter(j))); 
      sf_syst_tag_cuts->SetBinContent(nbin, TMath::Abs(
        eff_normal_tag->GetBinContent(nbin) / eff_mc_value - eff_alt_tag->GetBinContent(nbin) / eff_mc_alt_tag_value
      ));
      sf_syst_background_shape->SetBinContent(nbin, TMath::Abs(eff_syst_background_shape->GetBinContent(nbin) / eff_mc_value));
      sf_syst_signal_shape->SetBinContent(nbin, TMath::Abs(eff_syst_signal_shape->GetBinContent(nbin) / eff_mc_value));
      sf_syst_generator->SetBinContent(nbin, 
        TMath::Abs((eff_mc_value - eff_mc_value_LO)/(eff_mc_value)) * eff_data_value/eff_mc_value
      );
      //printf("Setting bin (%d,%d) of generator systematic to %f\n", i,j,syst_generator->GetBinContent(nbin));
      sf_syst_combined->SetBinContent(nbin, sqrt(
        pow(sf_syst_background_shape ->GetBinContent(nbin), 2)+
        pow(sf_syst_signal_shape     ->GetBinContent(nbin), 2)+
        pow(sf_syst_tag_cuts         ->GetBinContent(nbin), 2)+
        pow(sf_syst_generator        ->GetBinContent(nbin), 2)
      ));
    }}
  
    // make plots now
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.3f");
    mitPalette();
    TPaletteAxis *palette_axis;
  
    
    TCanvas *c_sf_syst_tag_cuts = new TCanvas("c_sf_syst_tag_cuts", "c_sf_syst_tag_cuts", 400, 600);
    sf_syst_tag_cuts->SetMinimum(0);
    sf_syst_tag_cuts->SetMaximum(0.2);
    sf_syst_tag_cuts->Draw("TEXT COLZ");
    sprintf(title, "Syst. unc. for %s. SF (%s WP) from choice of tag cuts", flavors_short, selections_caps[iwp]);
    sf_syst_tag_cuts->SetTitle(title);
    gPad->Update();
    sf_syst_tag_cuts->GetXaxis()->SetTitle("| #eta |");
    sf_syst_tag_cuts->GetXaxis()->SetTitleOffset(0.9);
    sf_syst_tag_cuts->GetXaxis()->SetTitleSize(0.04);
    sf_syst_tag_cuts->GetXaxis()->SetLabelSize(0.02);
    sf_syst_tag_cuts->GetYaxis()->SetTitle("p_{T} [GeV]");
    sf_syst_tag_cuts->GetYaxis()->SetTitleOffset(0.9);
    sf_syst_tag_cuts->GetYaxis()->SetTitleSize(0.04);
    sf_syst_tag_cuts->GetYaxis()->SetLabelSize(0.02);
    sf_syst_tag_cuts->GetYaxis()->SetRangeUser(10,100);
    sf_syst_tag_cuts->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) sf_syst_tag_cuts->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_sf_syst_tag_cuts->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_tag_cuts.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_sf_syst_tag_cuts->Print(filename);
    
    TCanvas *c_sf_syst_generator = new TCanvas("c_sf_syst_generator", "c_sf_syst_tag_cuts", 400, 600);
    sf_syst_generator->SetMinimum(0);
    sf_syst_generator->SetMaximum(0.2);
    sf_syst_generator->Draw("TEXT COLZ");
    gPad->Update();
    sprintf(title, "Syst. unc. for %s. SF (%s WP) from choice of generator", flavors_short, selections_caps[iwp]);
    sf_syst_generator->SetTitle(title);
    sf_syst_generator->GetXaxis()->SetTitle("| #eta |");
    sf_syst_generator->GetXaxis()->SetTitleOffset(0.9);
    sf_syst_generator->GetXaxis()->SetTitleSize(0.04);
    sf_syst_generator->GetXaxis()->SetLabelSize(0.02);
    sf_syst_generator->GetYaxis()->SetTitle("p_{T} [GeV]");
    sf_syst_generator->GetYaxis()->SetTitleOffset(0.9);
    sf_syst_generator->GetYaxis()->SetTitleSize(0.04);
    sf_syst_generator->GetYaxis()->SetLabelSize(0.02);
    sf_syst_generator->GetYaxis()->SetRangeUser(10,100);
    sf_syst_generator->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) sf_syst_generator->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_sf_syst_generator->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_generator.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_sf_syst_generator->Print(filename);
    
    TCanvas *c_sf_syst_signal_shape = new TCanvas("c_sf_syst_signal_shape", "c_sf_syst_signal_shape", 400, 600);
    sf_syst_signal_shape->SetMinimum(0);
    sf_syst_signal_shape->SetMaximum(0.2);
    sf_syst_signal_shape->Draw("TEXT COLZ");
    gPad->Update();
    sprintf(title, "Syst. unc. for %s. SF (%s WP) from signal shape", flavors_short, selections_caps[iwp]);
    sf_syst_signal_shape->SetTitle(title);
    sf_syst_signal_shape->GetXaxis()->SetTitle("| #eta |");
    sf_syst_signal_shape->GetXaxis()->SetTitleOffset(0.9);
    sf_syst_signal_shape->GetXaxis()->SetTitleSize(0.04);
    sf_syst_signal_shape->GetXaxis()->SetLabelSize(0.02);
    sf_syst_signal_shape->GetYaxis()->SetTitle("p_{T} [GeV]");
    sf_syst_signal_shape->GetYaxis()->SetTitleOffset(0.9);
    sf_syst_signal_shape->GetYaxis()->SetTitleSize(0.04);
    sf_syst_signal_shape->GetYaxis()->SetLabelSize(0.02);
    sf_syst_signal_shape->GetYaxis()->SetRangeUser(10,100);
    sf_syst_signal_shape->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) sf_syst_signal_shape->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_sf_syst_signal_shape->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_signal_shape.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_sf_syst_signal_shape->Print(filename);
    
    TCanvas *c_sf_syst_background_shape = new TCanvas("c_sf_syst_background_shape", "c_sf_syst_background_shape", 400, 600);
    sf_syst_background_shape->SetMinimum(0);
    sf_syst_background_shape->SetMaximum(0.2);
    sf_syst_background_shape->Draw("TEXT COLZ");
    gPad->Update();
    sprintf(title, "Syst. unc. for %s. SF (%s WP) from background shape", flavors_short, selections_caps[iwp]);
    sf_syst_background_shape->SetTitle(title);
    sf_syst_background_shape->GetXaxis()->SetTitle("| #eta |");
    sf_syst_background_shape->GetXaxis()->SetTitleOffset(0.9);
    sf_syst_background_shape->GetXaxis()->SetTitleSize(0.04);
    sf_syst_background_shape->GetXaxis()->SetLabelSize(0.02);
    sf_syst_background_shape->GetYaxis()->SetTitle("p_{T} [GeV]");
    sf_syst_background_shape->GetYaxis()->SetTitleOffset(0.9);
    sf_syst_background_shape->GetYaxis()->SetTitleSize(0.04);
    sf_syst_background_shape->GetYaxis()->SetLabelSize(0.02);
    sf_syst_background_shape->GetYaxis()->SetRangeUser(10,100);
    sf_syst_background_shape->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) sf_syst_background_shape->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_sf_syst_background_shape->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_background_shape.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_sf_syst_background_shape->Print(filename);
    
    TCanvas *c_sf_syst_combined = new TCanvas("c_sf_syst_combined", "c_sf_syst_combined", 400, 600);
    sf_syst_combined->SetMinimum(0);
    sf_syst_combined->SetMaximum(0.2);
    sf_syst_combined->Draw("TEXT COLZ");
    gPad->Update();
    sprintf(title, "Combined syst. unc. for %s. SF (%s WP)", flavors_short, selections_caps[iwp]);
    sf_syst_combined->SetTitle(title);
    sf_syst_combined->GetXaxis()->SetTitle("| #eta |");
    sf_syst_combined->GetXaxis()->SetTitleOffset(0.9);
    sf_syst_combined->GetXaxis()->SetTitleSize(0.04);
    sf_syst_combined->GetXaxis()->SetLabelSize(0.02);
    sf_syst_combined->GetYaxis()->SetTitle("p_{T} [GeV]");
    sf_syst_combined->GetYaxis()->SetTitleOffset(0.9);
    sf_syst_combined->GetYaxis()->SetTitleSize(0.04);
    sf_syst_combined->GetYaxis()->SetLabelSize(0.02);
    sf_syst_combined->GetYaxis()->SetRangeUser(10,100);
    sf_syst_combined->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) sf_syst_combined->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_sf_syst_combined->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_combined.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_sf_syst_combined->Print(filename);
    // save to file
    sf_file->cd();
    sf_file->Delete((string(sf_syst_background_shape ->GetName()) + ";*").c_str()); 
    sf_file->Delete((string(sf_syst_signal_shape     ->GetName()) + ";*").c_str()); 
    sf_file->Delete((string(sf_syst_generator        ->GetName()) + ";*").c_str()); 
    sf_file->Delete((string(sf_syst_tag_cuts         ->GetName()) + ";*").c_str());   
    sf_file->Delete((string(sf_syst_combined         ->GetName()) + ";*").c_str());     
    sf_syst_background_shape ->Write(); 
    sf_syst_signal_shape     ->Write(); 
    sf_syst_generator        ->Write(); 
    sf_syst_tag_cuts         ->Write();   
    sf_syst_combined         ->Write();     

  } 
  sf_file->Close();
}

void check_toy() {
  printf("wow\n");
}
/*void combined_toy_systematics_muons(string sf_filename, string plots_dir, string root_dir) {
  unsigned int n_working_points=3;
  char const *flavors_short[]        = {"mu"        , "mu"        , "mu"       };
  char const *flavors_caps[]         = {"Muon"      , "Muon"      , "Muon"     };
  char const *flavors_lowercase[]    = {"muon"      , "muon"      , "muon"     };
  char const *selections_caps[]      = {"Loose"     , "Medium"    , "Tight"    };
  char const *selections_lowercase[] = {"loose"     , "medium"    , "tight"    };
  char const directory_format[]="Single%s_BaselineTo%s%s_%sTnP";
  char const diff_tag_str[]="_diffTag";

  if( plots_dir[plots_dir.size()-1]  != '/' ) plots_dir = plots_dir + "/";
  if( root_dir[root_dir.size()-1]  != '/' ) root_dir = root_dir + "/";
  TFile *sf_file = TFile::Open(sf_filename.c_str(), "UPDATE");
  TFile *tag_cuts_file, *background_file, *signal_file, *generator_file;
  TH2D *eff_mc, *eff_mc_LO, *eff_mc_diffTag, *eff_data;
  TFile *file_eff_mc_LO, *file_eff_mc_diffTag;
  char filename[512], object_name[128], title[512];
  for(unsigned int iwp=0; iwp<n_working_points; iwp++) {
    sprintf(filename, "~/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/DYJetsToLL_LO_BaselineTo%s_muonTnP/eff.root", selections_caps[iwp]);
    file_eff_mc_LO = TFile::Open(filename);
    eff_mc_LO = (TH2D*) file_eff_mc_LO->Get("hEffEtaPt");
    sprintf(filename, "~/TagAndProbe/2016-05-02_76x_muons_coarseEta/template_erfcexp/DYJetsToLL_BaselineTo%s_diffTag_%sTnP/eff.root", selections_caps[iwp], flavors_lowercase);
    file_eff_mc_diffTag=TFile::Open(filename);
    eff_mc_diffTag = (TH2D*)file_eff_mc_diffTag->Get("hEffEtaPt");

    sprintf(object_name, "eff_mc_%s_%s", selections_caps[iwp], flavors_short); 
    eff_mc = (TH2D*) sf_file->Get(object_name);
    sprintf(object_name, "eff_data_%s_%s", selections_caps[iwp], flavors_short);
    eff_data = (TH2D*) sf_file->Get(object_name);
  
    sprintf(filename, "%sSingle%s_BaselineTo%s_diffTag_%sTnP/toystudy_template_erfcexp_diffTag_versus_template_erfcexp_diffTag.root", root_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    tag_cuts_file   = TFile::Open(filename);
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/toystudy_template_erfcexp_versus_template_exp.root", root_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    background_file = TFile::Open(filename);
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/toystudy_template_erfcexp_versus_BWCBPlusVoigt_erfcexp.root", root_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    signal_file     = TFile::Open(filename);
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/toystudy_template_erfcexp_versus_template_erfcexp_LO.root", root_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    generator_file  = TFile::Open(filename);
  
    sprintf(object_name, "scalefactors_%s_%s_syst_error_background_shape", selections_caps[iwp], flavors_caps);
    TH2D *syst_background_shape = (TH2D*) ((TH2D*) background_file->Get("syst_from_Ztest"))->Clone(object_name);
    sprintf(object_name, "scalefactors_%s_%s_syst_error_signal_shape", selections_caps[iwp], flavors_caps);
    TH2D *syst_signal_shape     = (TH2D*) ((TH2D*) signal_file->Get("syst_from_Ztest"))->Clone(object_name);
    sprintf(object_name, "scalefactors_%s_%s_syst_generator", selections_caps[iwp], flavors_caps);
    TH2D *syst_generator        = (TH2D*) ((TH2D*) generator_file->Get("syst_from_Ztest"))->Clone(object_name);

    TH2D *eff_normal_tag = (TH2D*) signal_file->Get("mean_eff_method1");
    TH2D *eff_alt_tag    = (TH2D*) tag_cuts_file->Get("mean_eff_method1");
    sprintf(object_name, "scalefactors_%s_%s_syst_error_tag_cuts", selections_caps[iwp], flavors_caps);
    TH2D *syst_tag_cuts  = (TH2D*) eff_normal_tag->Clone(object_name);
    syst_tag_cuts->Add(eff_alt_tag, -1);
    
    sprintf(object_name, "scalefactors_%s_%s_syst_error_combined", selections_caps[iwp], flavors_caps);
    TH2D *syst_combined = (TH2D*) syst_generator->Clone(object_name);
  
    for(unsigned int i = 1; i <= (unsigned int) syst_tag_cuts->GetNbinsX(); i++) { for(unsigned int j = 1; j <= (unsigned int) syst_tag_cuts->GetNbinsY(); j++) {
      int nbin = syst_tag_cuts->GetBin(i,j);
      double eff_data_value    = eff_data->GetBinContent(eff_data->FindBin(syst_tag_cuts->GetXaxis()->GetBinCenter(i), syst_tag_cuts->GetYaxis()->GetBinCenter(j)));
      double eff_mc_value    = eff_mc->GetBinContent(eff_mc->FindBin(syst_tag_cuts->GetXaxis()->GetBinCenter(i), syst_tag_cuts->GetYaxis()->GetBinCenter(j)));
      double eff_mc_diffTag_value    = eff_mc_diffTag->GetBinContent(eff_mc_diffTag->FindBin(syst_tag_cuts->GetXaxis()->GetBinCenter(i), syst_tag_cuts->GetYaxis()->GetBinCenter(j)));
      double eff_mc_value_LO = eff_mc_LO->GetBinContent(eff_mc_LO->FindBin(syst_tag_cuts->GetXaxis()->GetBinCenter(i), syst_tag_cuts->GetYaxis()->GetBinCenter(j))); 
      syst_tag_cuts->SetBinContent(nbin, TMath::Abs(
        eff_normal_tag->GetBinContent(nbin) / eff_mc_value - eff_alt_tag->GetBinContent(nbin) / eff_mc_diffTag_value
      ));
      syst_background_shape->SetBinContent(nbin, syst_background_shape->GetBinContent(nbin) / eff_mc_value);
      syst_signal_shape->SetBinContent(nbin, syst_signal_shape->GetBinContent(nbin) / eff_mc_value);
      syst_generator->SetBinContent(nbin, TMath::Abs(eff_data_value / eff_mc_value - eff_data_value / eff_mc_value_LO));
      printf("Setting bin (%d,%d) of generator systematic to %f\n", i,j,syst_generator->GetBinContent(nbin));
      syst_combined->SetBinContent(nbin, sqrt(
        pow(syst_background_shape ->GetBinContent(nbin), 2)+
        pow(syst_signal_shape     ->GetBinContent(nbin), 2)+
        pow(syst_tag_cuts         ->GetBinContent(nbin), 2)+
        pow(syst_generator        ->GetBinContent(nbin), 2)
      ));
    }}
  
    // make plots now
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.3f");
    mitPalette();
    TPaletteAxis *palette_axis;
  
    
    TCanvas *c_syst_tag_cuts = new TCanvas("c_syst_tag_cuts", "c_syst_tag_cuts", 400, 600);
    sprintf(title, "Syst. unc. for %s. SF (%s WP) from choice of tag cuts", flavors_short, selections_caps[iwp]);
    syst_tag_cuts->SetTitle(title);
    syst_tag_cuts->Draw("TEXT COLZ");
    gPad->Update();
    syst_tag_cuts->GetXaxis()->SetTitle("| #eta |");
    syst_tag_cuts->GetXaxis()->SetTitleOffset(0.9);
    syst_tag_cuts->GetXaxis()->SetTitleSize(0.04);
    syst_tag_cuts->GetXaxis()->SetLabelSize(0.02);
    syst_tag_cuts->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_tag_cuts->GetYaxis()->SetTitleOffset(0.9);
    syst_tag_cuts->GetYaxis()->SetTitleSize(0.04);
    syst_tag_cuts->GetYaxis()->SetLabelSize(0.02);
    syst_tag_cuts->GetYaxis()->SetRangeUser(10,100);
    syst_tag_cuts->SetMinimum(0);
    syst_tag_cuts->SetMaximum(0.2);
    syst_tag_cuts->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_tag_cuts->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_tag_cuts->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_tag_cuts.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_syst_tag_cuts->Print(filename);
    
    TCanvas *c_syst_generator = new TCanvas("c_syst_generator", "c_syst_tag_cuts", 400, 600);
    syst_generator->Draw("TEXT COLZ");
    gPad->Update();
    sprintf(title, "Syst. unc. for %s. SF (%s WP) from choice of generator", flavors_short, selections_caps[iwp]);
    syst_generator->SetTitle(title);
    syst_generator->GetXaxis()->SetTitle("| #eta |");
    syst_generator->GetXaxis()->SetTitleOffset(0.9);
    syst_generator->GetXaxis()->SetTitleSize(0.04);
    syst_generator->GetXaxis()->SetLabelSize(0.02);
    syst_generator->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_generator->GetYaxis()->SetTitleOffset(0.9);
    syst_generator->GetYaxis()->SetTitleSize(0.04);
    syst_generator->GetYaxis()->SetLabelSize(0.02);
    syst_generator->GetYaxis()->SetRangeUser(10,100);
    syst_generator->SetMinimum(0);
    syst_generator->SetMaximum(0.2);
    syst_generator->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_generator->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_generator->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_generator.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_syst_generator->Print(filename);
    
    TCanvas *c_syst_signal_shape = new TCanvas("c_syst_signal_shape", "c_syst_signal_shape", 400, 600);
    syst_signal_shape->Draw("TEXT COLZ");
    gPad->Update();
    sprintf(title, "Syst. unc. for %s. SF (%s WP) from signal shape", flavors_short, selections_caps[iwp]);
    syst_signal_shape->SetTitle(title);
    syst_signal_shape->GetXaxis()->SetTitle("| #eta |");
    syst_signal_shape->GetXaxis()->SetTitleOffset(0.9);
    syst_signal_shape->GetXaxis()->SetTitleSize(0.04);
    syst_signal_shape->GetXaxis()->SetLabelSize(0.02);
    syst_signal_shape->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_signal_shape->GetYaxis()->SetTitleOffset(0.9);
    syst_signal_shape->GetYaxis()->SetTitleSize(0.04);
    syst_signal_shape->GetYaxis()->SetLabelSize(0.02);
    syst_signal_shape->GetYaxis()->SetRangeUser(10,100);
    syst_signal_shape->SetMinimum(0);
    syst_signal_shape->SetMaximum(0.2);
    syst_signal_shape->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_signal_shape->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_signal_shape->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_signal_shape.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_syst_signal_shape->Print(filename);
    
    TCanvas *c_syst_background_shape = new TCanvas("c_syst_background_shape", "c_syst_background_shape", 400, 600);
    syst_background_shape->Draw("TEXT COLZ");
    gPad->Update();
    sprintf(title, "Syst. unc. for %s. SF (%s WP) from background shape", flavors_short, selections_caps[iwp]);
    syst_background_shape->SetTitle(title);
    syst_background_shape->GetXaxis()->SetTitle("| #eta |");
    syst_background_shape->GetXaxis()->SetTitleOffset(0.9);
    syst_background_shape->GetXaxis()->SetTitleSize(0.04);
    syst_background_shape->GetXaxis()->SetLabelSize(0.02);
    syst_background_shape->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_background_shape->GetYaxis()->SetTitleOffset(0.9);
    syst_background_shape->GetYaxis()->SetTitleSize(0.04);
    syst_background_shape->GetYaxis()->SetLabelSize(0.02);
    syst_background_shape->GetYaxis()->SetRangeUser(10,100);
    syst_background_shape->SetMinimum(0);
    syst_background_shape->SetMaximum(0.2);
    syst_background_shape->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_background_shape->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_background_shape->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_background_shape.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_syst_background_shape->Print(filename);
    
    TCanvas *c_syst_combined = new TCanvas("c_syst_combined", "c_syst_combined", 400, 600);
    syst_combined->Draw("TEXT COLZ");
    gPad->Update();
    sprintf(title, "Combined syst. unc. for %s. SF (%s WP)", flavors_short, selections_caps[iwp]);
    syst_combined->SetTitle(title);
    syst_combined->GetXaxis()->SetTitle("| #eta |");
    syst_combined->GetXaxis()->SetTitleOffset(0.9);
    syst_combined->GetXaxis()->SetTitleSize(0.04);
    syst_combined->GetXaxis()->SetLabelSize(0.02);
    syst_combined->GetYaxis()->SetTitle("p_{T} [GeV]");
    syst_combined->GetYaxis()->SetTitleOffset(0.9);
    syst_combined->GetYaxis()->SetTitleSize(0.04);
    syst_combined->GetYaxis()->SetLabelSize(0.02);
    syst_combined->GetYaxis()->SetRangeUser(10,100);
    syst_combined->SetMinimum(0);
    syst_combined->SetMaximum(0.2);
    syst_combined->SetMarkerSize(1.3);
    palette_axis = (TPaletteAxis*) syst_combined->GetListOfFunctions()->FindObject("palette"); 
    palette_axis->SetLabelSize(0.02);
    c_syst_combined->Update();
    sprintf(filename, "%sSingle%s_BaselineTo%s_%sTnP/syst_combined.png", plots_dir.c_str(), flavors_caps, selections_caps[iwp], flavors_lowercase);
    c_syst_combined->Print(filename);
    // save to file
    sf_file->cd();
    syst_background_shape ->Write(); 
    syst_signal_shape     ->Write(); 
    syst_generator        ->Write(); 
    syst_tag_cuts         ->Write();   
    syst_combined         ->Write();     

  } 
  sf_file->Close();
}*/

void do_the_electrons() {

  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_diffTag_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_diffTag_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_diffTag_electronTnP/", "template_erfcexp_diffTag", "template_erfcexp_diffTag", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "template_erfcexp", "template_exp", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "template_erfcexp", "template_erfcexp_LO", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToVeto_electronTnP/", "template_erfcexp", "BWCBPlusVoigt_erfcexp", 1000);

  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_diffTag_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_diffTag_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_diffTag_electronTnP/", "template_erfcexp_diffTag", "template_erfcexp_diffTag", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "template_erfcexp", "template_exp", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "template_erfcexp", "template_erfcexp_LO", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToLoose_electronTnP/", "template_erfcexp", "BWCBPlusVoigt_erfcexp", 1000);

  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_diffTag_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_diffTag_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_diffTag_electronTnP/", "template_erfcexp_diffTag", "template_erfcexp_diffTag", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "template_erfcexp", "template_exp", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "template_erfcexp", "template_erfcexp_LO", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToMedium_electronTnP/", "template_erfcexp", "BWCBPlusVoigt_erfcexp", 1000);

  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_diffTag_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_diffTag_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_diffTag_electronTnP/", "template_erfcexp_diffTag", "template_erfcexp_diffTag", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "template_erfcexp", "template_exp", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "template_erfcexp", "template_erfcexp_LO", 1000);
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "plots/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "root/2016-05-02/2016-05-02_76x_electrons/template_erfcexp/SingleElectron_BaselineToTight_electronTnP/", "template_erfcexp", "BWCBPlusVoigt_erfcexp", 1000);
}
void do_the_muons() {
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_diffTag_muonTnP/", "plots/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_diffTag_muonTnP/", "root/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_diffTag_muonTnP/", "template_erfcexp_diffTag", "template_erfcexp_diffTag", 1000, "toy%06d", "muons");
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_muonTnP/", "plots/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_muonTnP/", "root/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_muonTnP/", "template_erfcexp", "template_exp", 1000, "toy%06d", "muons");
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_muonTnP/", "plots/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_muonTnP/", "root/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToLoose_muonTnP/", "template_erfcexp", "BWCBPlusVoigt_erfcexp", 1000, "toy%06d", "muons");

  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_diffTag_muonTnP/", "plots/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_diffTag_muonTnP/", "root/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_diffTag_muonTnP/", "template_erfcexp_diffTag", "template_erfcexp_diffTag", 1000, "toy%06d", "muons");
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_muonTnP/", "plots/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_muonTnP/", "root/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_muonTnP/", "template_erfcexp", "template_exp", 1000, "toy%06d", "muons");
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_muonTnP/", "plots/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_muonTnP/", "root/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToMedium_muonTnP/", "template_erfcexp", "BWCBPlusVoigt_erfcexp", 1000, "toy%06d", "muons");

  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_diffTag_muonTnP/", "plots/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_diffTag_muonTnP/", "root/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_diffTag_muonTnP/", "template_erfcexp_diffTag", "template_erfcexp_diffTag", 1000, "toy%06d", "muons");
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_muonTnP/", "plots/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_muonTnP/", "root/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_muonTnP/", "template_erfcexp", "template_exp", 1000, "toy%06d", "muons");
  compare_methods("/home/dhsu/TagAndProbe/toys/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_muonTnP/", "plots/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_muonTnP/", "root/2016-05-02/2016-05-02_76x_muons_coarseEta/template_erfcexp/SingleMuon_BaselineToTight_muonTnP/", "template_erfcexp", "BWCBPlusVoigt_erfcexp", 1000, "toy%06d", "muons");
}
