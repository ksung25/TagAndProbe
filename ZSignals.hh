#include "TROOT.h"
#include "TH1D.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooVoigtianShape.h"
#include "RooVoigtian.h"
#include "RooKeysPdf.h"


class CSignalModel
{
public:
  CSignalModel():model(0){}
  virtual ~CSignalModel(){ delete model; }
  RooAbsPdf *model;
protected:
  std::vector<double> readSigParams(
    const int ibin,
    const std::string name,
    std::vector<std::string> paramNames,
    std::string refDir
  );

};

class CBreitWignerConvCrystalBall : public CSignalModel
{
public:
  CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass);
  ~CBreitWignerConvCrystalBall();
  RooRealVar     *mass, *width;
  RooBreitWigner *bw;
  RooRealVar     *mean, *sigma, *alpha, *n;
  RooCBShape     *cb;
};

class CBWCBPlusVoigt : public CSignalModel
{
public:
  CBWCBPlusVoigt(RooRealVar &m, const Bool_t pass, double ptMin, double ptMax, int ibin=0, const std::string fitname="", std::string refDir="");
  ~CBWCBPlusVoigt();
  RooRealVar     *mass, *width;
  RooBreitWigner *bw;
  RooRealVar     *mean, *sigma, *alpha, *n, *vMean, *vWidth, *vSigma, *fsrFrac;
  RooFormulaVar  *oneMinusFsrFrac;
  RooCBShape     *cb;
  RooVoigtian    *voigt;
  RooAbsPdf      *bwcb;
};

class CMCTemplateConvGaussian : public CSignalModel
{
public:
  CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, RooRealVar *sigma0=0, int intOrder=1);
  ~CMCTemplateConvGaussian();
  RooRealVar  *mean, *sigma;
  RooGaussian *gaus;
  TH1D        *inHist;
  RooDataHist *dataHist;
  RooHistPdf  *histPdf;
};

class CVoigtianCBShape : public CSignalModel
{
public:
  CVoigtianCBShape(RooRealVar &m, const Bool_t pass);
  ~CVoigtianCBShape();
  RooRealVar *mass, *width;
  RooRealVar *sigma, *alpha, *n;
};

class CMCDatasetConvGaussian : public CSignalModel
{
public:
  CMCDatasetConvGaussian(RooRealVar &m, TTree* tree, const Bool_t pass, RooRealVar *sigma0=0);
  ~CMCDatasetConvGaussian();
  RooRealVar  *mean, *sigma;
  RooGaussian *gaus;
  TTree       *inTree;
  RooDataSet  *dataSet;
  RooKeysPdf  *keysPdf;
};

//--------------------------------------------------------------------------------------------------
CBreitWignerConvCrystalBall::CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  
  sprintf(vname,"mass%s",name);
  mass = new RooRealVar(vname,vname,91,80,100);    
  mass->setVal(91.1876);
  mass->setConstant(kTRUE);
  
  sprintf(vname,"width%s",name);
  width = new RooRealVar(vname,vname,2.5,0.1,10);    
  width->setVal(2.4952);
  width->setConstant(kTRUE);
  
  sprintf(vname,"bw%s",name);
  bw = new RooBreitWigner(vname,vname,m,*mass,*width);

  if(pass) {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,5);
    sprintf(vname,"alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
    sprintf(vname,"n%s",name);     n     = new RooRealVar(vname,vname,1,0,10);
  } else {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,5);
    sprintf(vname,"alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
    sprintf(vname,"n%s",name);     n     = new RooRealVar(vname,vname,1,0,10);
  }  
//  n->setVal(1.0);
//  n->setConstant(kTRUE);
  
  sprintf(vname,"cb%s",name);
  cb = new RooCBShape(vname,vname,m,*mean,*sigma,*alpha,*n);
        
  sprintf(vname,"signal%s",name);
  model = new RooFFTConvPdf(vname,vname,m,*bw,*cb);
}

CBreitWignerConvCrystalBall::~CBreitWignerConvCrystalBall()
{
  delete mass;  mass=0;
  delete width; width=0;
  delete bw;    bw=0;
  delete mean;  mean=0;
  delete sigma; sigma=0;
  delete alpha; alpha=0;
  delete n;     n=0;
  delete cb;    cb=0;
}

//--------------------------------------------------------------------------------------------------
CBWCBPlusVoigt::CBWCBPlusVoigt(RooRealVar &m, const Bool_t pass, double ptMin, double ptMax, const int ibin, const std::string fitname, std::string refDir)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  char formula[100];
  
  sprintf(vname,"mass%s",name);
  mass = new RooRealVar(vname,vname,91,80,100);    
  mass->setVal(91.1876);
  mass->setConstant(kTRUE);
  
  sprintf(vname,"width%s",name);
  width = new RooRealVar(vname,vname,2.5,0.1,10);    
  width->setVal(2.4952);
  width->setConstant(kTRUE);
  
  sprintf(vname,"bw%s",name);
  bw = new RooBreitWigner(vname,vname,m,*mass,*width);
  
  double fsrPeak=80;
  double fsrSigma=7;
  double fsrSigmaMax=10;
  double fsrSigmaMin=6;
  double fsrPeakRange=20;
  double fsrFracMin=0.05;
  double fsrFracMax=0.2;
  double fsrFracInit=0.05;
  //fsrFracInit = 0.2;
  //fsrFracMax = 0.3;
  if(ptMin >= 100) {
    fsrFracMin=0;
    fsrPeak=140;
  //} else if(ptMin >= 70) {
  //  fsrPeak=115;
  } else if(ptMin >= 50) {
    fsrPeak=105;
    fsrPeakRange=10;
    //fsrFracInit=0.2;
    fsrFracMin=0;
    //fsrFracMax = 0.3;
    //fsrSigma=7;
  } else if(ptMin >= 40) {
    fsrPeak=95;
    //fsrFracMin=.1;
    //fsrFracMax = 0.3;
  } else if(ptMin >= 30) {
    fsrPeak=80;
    //fsrFracMax = 0.3;
    //fsrFracMin=.1;
    //fsrPeakRange=10;
  } else if(ptMin >= 20) {
    fsrPeak=75;
    //fsrFracMin=.15;
  } else if(ptMin >= 10) {
    fsrPeak=55;
    //fsrFracMin=.3;
    //fsrFracInit=.4;
    //fsrFracMax=.8;
    fsrPeakRange=5;
    //fsrSigma=10;
  } else {
    fsrPeakRange=300;
  }

  //sprintf(vname,"vMean%s",name);      vMean  = new RooRealVar(vname,vname,60,5,150);
  //sprintf(vname,"vMean%s",name);      vMean  = new RooRealVar(vname,vname,fsrPeak,5,150);
  sprintf(vname,"vMean%s",name);      vMean  = new RooRealVar(vname,vname,fsrPeak,fsrPeak-fsrPeakRange,fsrPeak+fsrPeakRange);
  sprintf(vname,"vWidth%s",name);     vWidth = new RooRealVar(vname,vname,0.01,0.001,.1);
  sprintf(vname,"vSigma%s",name);     vSigma = new RooRealVar(vname,vname,fsrSigma, fsrSigmaMin, fsrSigmaMax);
  if(pass) {
    sprintf(vname,"mean%s",name);       mean   = new RooRealVar(vname,vname,0,-5,5);
    sprintf(vname,"sigma%s",name);      sigma  = new RooRealVar(vname,vname,2,1,5);
    sprintf(vname,"alpha%s",name);      alpha  = new RooRealVar(vname,vname,.8, -1, 3);
    sprintf(vname,"n%s",name);          n      = new RooRealVar(vname,vname,1.5,0,10);
    sprintf(vname,"fsrFrac%s",name);    fsrFrac= new RooRealVar(vname,vname, 0.01, 0,fsrFracMax);
  } else {
    sprintf(vname,"mean%s",name);       mean  = new RooRealVar(vname,vname,0,-5,5);
    sprintf(vname,"sigma%s",name);      sigma = new RooRealVar(vname,vname,2,1,5);
    sprintf(vname,"alpha%s",name);      alpha = new RooRealVar(vname,vname,.8,-1,3);
    sprintf(vname,"n%s",name);          n     = new RooRealVar(vname,vname,1.5,0.5,10);
    sprintf(vname,"fsrFrac%s",name);    fsrFrac= new RooRealVar(vname,vname, fsrFracInit, fsrFracMin,fsrFracMax);
  }
  sprintf(formula, "1 - fsrFrac%s", name);
  sprintf(vname,"oneMinusFsrFrac%s",name);    oneMinusFsrFrac= new RooFormulaVar(vname,vname,formula, *fsrFrac);
  if(refDir != "") {
    std::vector<std::string> paramNames = {
      "vMean",
      "vWidth",
      "vSigma",
      "mean",
      "sigma",
      "alpha",
      "n",
      "fsrFrac"
    };
    for(unsigned int i=0; i<paramNames.size(); i++) {
      if(pass) {
        paramNames[i] = paramNames[i] + "Pass";
      } else {
        paramNames[i] = paramNames[i] + "Fail";
      }
    }
    std::vector<double> params = readSigParams(ibin, fitname, paramNames, refDir);
    vMean    ->setVal(params[0]);  vMean    ->setMin(.9*params[0]);  vMean    ->setMax(1.1*params[0]);  
    vWidth   ->setVal(params[1]);    
    vSigma   ->setVal(params[2]);  vSigma   ->setMin(.9*params[2]);  vSigma   ->setMax(1.1*params[2]);  
    mean     ->setVal(params[3]);  
    sigma    ->setVal(params[4]);  
    alpha    ->setVal(params[5]);  alpha    ->setMin(params[5]-2);   alpha    ->setMax(params[5]+2);
    n        ->setVal(TMath::Min(10., TMath::Max(0., params[6])));  
    n        ->setMin(TMath::Max(0., .9*params[6]));  
    n        ->setMax(TMath::Min(10., 1.1*params[6]));
    fsrFrac  ->setVal(params[7]);  fsrFrac  ->setMin(.9*params[7]);  fsrFrac  ->setMax(1.1*params[7]);  

  }
  
   
  sprintf(vname,"cb%s",name);
  cb = new RooCBShape(vname,vname,m,*mean,*sigma,*alpha,*n);
  
  sprintf(vname, "voigt%s",name);
  voigt = new RooVoigtian(vname,vname, m, *vMean, *vWidth, *vSigma);
  sprintf(vname,"bwcb%s",name);
  bwcb = new RooFFTConvPdf(vname,vname,m,*bw,*cb);
  RooArgList *pdfs = new RooArgList(*voigt, *bwcb);
  RooArgList *coeffs = new RooArgList(*fsrFrac, *oneMinusFsrFrac); 
  sprintf(vname,"signal%s",name);
  model = new RooAddPdf(vname, vname, *pdfs, *coeffs);
  //model = new RooFFTConvPdf(vname,vname,m,*bw,*cb);
}

CBWCBPlusVoigt::~CBWCBPlusVoigt()
{
  delete vMean;    vMean=0;
  delete vWidth;    vWidth=0;
  delete vSigma;    vSigma=0;
  delete oneMinusFsrFrac;    oneMinusFsrFrac=0;
  delete fsrFrac;    fsrFrac=0;
  delete voigt; voigt=0;
  delete bwcb;  bwcb=0;
  delete mass;  mass=0;
  delete width; width=0;
  delete bw;    bw=0;
  delete mean;  mean=0;
  delete sigma; sigma=0;
  delete alpha; alpha=0;
  delete n;     n=0;
  delete cb;    cb=0;
}

//--------------------------------------------------------------------------------------------------
CMCTemplateConvGaussian::CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, RooRealVar *sigma0, int intOrder)
{  
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];  
  
  if(pass) {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,-1,-5,5);
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,0.5,0.5,5); }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  } else {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,-1,-5,5);
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,0.5,0.5,5); }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  }


  sprintf(vname,"inHist_%s",hist->GetName());
  inHist = (TH1D*)hist->Clone(vname);
  
  sprintf(vname,"dataHist%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(m),inHist);
  sprintf(vname,"histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,intOrder);
  sprintf(vname,"signal%s",name);   model    = new RooFFTConvPdf(vname,vname,m,*histPdf,*gaus);
}

CMCTemplateConvGaussian::~CMCTemplateConvGaussian()
{
  delete mean;     mean=0;
  //delete sigma;    sigma=0;
  delete gaus;     gaus=0;
  delete inHist;   inHist=0;
  delete dataHist; dataHist=0;
  delete histPdf;  histPdf=0;
}

//--------------------------------------------------------------------------------------------------
CVoigtianCBShape::CVoigtianCBShape(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  
  sprintf(vname,"mass%s",name);
  mass = new RooRealVar(vname,vname,91,80,100);
    
  sprintf(vname,"width%s",name);
  width = new RooRealVar(vname,vname,2.5,0.1,10);    
  width->setVal(2.4952);
  width->setConstant(kTRUE);
  
  if(pass) {
    sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,10);
    sprintf(vname,"alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
  } else {
    sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,10);
    sprintf(vname,"alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
  }
  
  sprintf(vname,"n%s",name);     
  n = new RooRealVar(vname,vname,1,0,10);
  n->setVal(1.0);
  
  sprintf(vname,"signal%s",name);
  model = new RooVoigtianShape(vname,vname,m,*mass,*sigma,*alpha,*n,*width,0);  
}

CVoigtianCBShape::~CVoigtianCBShape()
{
  delete mass;  mass=0;
  delete width; width=0;
  delete sigma; sigma=0;
  delete alpha; alpha=0;
  delete n;     n=0;
}

//--------------------------------------------------------------------------------------------------
CMCDatasetConvGaussian::CMCDatasetConvGaussian(RooRealVar &m, TTree* tree, const Bool_t pass, RooRealVar *sigma0)
{  
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];  

  if(pass) {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,5); }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  } else {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,5); }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  }
  
  sprintf(vname,"inTree_%s",tree->GetName());
  inTree = (TTree*)tree->Clone(vname);
  
  sprintf(vname,"dataSet%s",name); dataSet = new RooDataSet(vname,vname,inTree,RooArgSet(m));
  sprintf(vname,"keysPdf%s",name); keysPdf = new RooKeysPdf(vname,vname,m,*dataSet,RooKeysPdf::NoMirror,1);
  sprintf(vname,"signal%s",name);  model   = new RooFFTConvPdf(vname,vname,m,*keysPdf,*gaus);
}

CMCDatasetConvGaussian::~CMCDatasetConvGaussian()
{
  delete mean;    mean=0;
  delete sigma;   sigma=0;
  delete gaus;    gaus=0;
  delete inTree;  inTree=0;
  delete dataSet; dataSet=0;
  delete keysPdf; keysPdf=0;
}
//--------------------------------------------------------------------------------------------------
std::vector<double> CSignalModel::readSigParams(
  const int ibin,
  const std::string fitname,
  std::vector<std::string> paramNames,
  std::string refDir
) {
  std::vector<double> params;
  char rname[512];
  sprintf(
    rname,
    "%s/plots/fitres%s_%i.txt",
    refDir.c_str(),
    fitname.c_str(),
    ibin
  );
  ifstream rfile;
  for(unsigned int i = 0; i<paramNames.size(); i++) {
    rfile.open(rname);
    assert(rfile.is_open());
    std::string line;
    bool found_param = false;
    printf("Looking for parameter %s...\n", paramNames[i].c_str());
    while(getline(rfile,line)) {
      //printf("%s\n", line.c_str());
      size_t found = line.find(" "+paramNames[i]+" ");
      if(found!=string::npos) {
        std::string varname, initval, finalval, pmstr, error, corr;
        std::stringstream ss(line);
        ss >> varname >> initval >> finalval >> pmstr >> error >> corr;
        params.push_back(atof(finalval.c_str()));
        printf("Got value %f\n", atof(finalval.c_str()));
        found_param=true;
        break;
      }
    }
    assert(found_param);
    rfile.close();
  }
  return params;
}
