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
  CBWCBPlusVoigt(RooRealVar &m, const Bool_t pass, double fsrPeak);
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
CBWCBPlusVoigt::CBWCBPlusVoigt(RooRealVar &m, const Bool_t pass, double fsrPeak)
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

  //sprintf(vname,"vMean%s",name);      vMean  = new RooRealVar(vname,vname,60,5,150);
  sprintf(vname,"vMean%s",name);      vMean  = new RooRealVar(vname,vname,fsrPeak,5,150);
  sprintf(vname,"vWidth%s",name);     vWidth = new RooRealVar(vname,vname,5,0.1,10);
  sprintf(vname,"vSigma%s",name);     vSigma = new RooRealVar(vname,vname,10,1,50);
  if(pass) {
    sprintf(vname,"mean%s",name);       mean   = new RooRealVar(vname,vname,0,-10,10);
    sprintf(vname,"sigma%s",name);      sigma  = new RooRealVar(vname,vname,1,0.1,15);
    sprintf(vname,"alpha%s",name);      alpha  = new RooRealVar(vname,vname,5,-20,20);
    sprintf(vname,"n%s",name);          n      = new RooRealVar(vname,vname,1,0,10);
    sprintf(vname,"fsrFrac%s",name);    fsrFrac= new RooRealVar(vname,vname,0.3,0,.45);
  } else {
    sprintf(vname,"mean%s",name);       mean  = new RooRealVar(vname,vname,0,-10,10);
    sprintf(vname,"sigma%s",name);      sigma = new RooRealVar(vname,vname,1,0.1,15);
    sprintf(vname,"alpha%s",name);      alpha = new RooRealVar(vname,vname,5,-20,20);
    sprintf(vname,"n%s",name);          n     = new RooRealVar(vname,vname,1,0,10);
    sprintf(vname,"fsrFrac%s",name);    fsrFrac= new RooRealVar(vname,vname,0.3,0,.45);
  }
  sprintf(formula, "1 - fsrFrac%s", name);
  sprintf(vname,"oneMinusFsrFrac%s",name);    oneMinusFsrFrac= new RooFormulaVar(vname,vname,formula, *fsrFrac);
  //sprintf(vname, "fsrNorm%s", name);    fsrNorm = new RooFormulaVar(vname, "fsrFrac", fsrFrac);
//  n->setVal(1.0);
//  n->setConstant(kTRUE);
  
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
