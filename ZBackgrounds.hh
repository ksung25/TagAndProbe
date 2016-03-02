#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooCMSShape.h"
#include "RooGenericPdf.h"
#include <vector>

class CBackgroundModel
{
public:
  CBackgroundModel():model(0){}
  virtual ~CBackgroundModel() { delete model; }
  RooAbsPdf *model;
protected:
  std::vector<double> readBkgParams(
    const int ibin,
    const std::string name,
    std::vector<std::string> paramNames,
    std::string refDir
  );
};

class CExponential : public CBackgroundModel
{
public:
  CExponential(RooRealVar &m, const Bool_t pass);
  ~CExponential();
  RooRealVar *t;
  RooRealVar *offset;
  RooFormulaVar *mMinusOffset;
};

class CErfcExpo : public CBackgroundModel
{
public:
  CErfcExpo(RooRealVar &m, const Bool_t pass);
  ~CErfcExpo();
  RooRealVar *alfa, *beta, *gamma, *peak; 
};

class CErfcExpoFixed : public CBackgroundModel
{
public:
  CErfcExpoFixed(RooRealVar &m, const Bool_t pass, const int ibin, const std::string name, std::string refDir );
  ~CErfcExpoFixed();
  RooRealVar *alfa, *beta, *gamma, *peak; 
};

class CDoubleExp : public CBackgroundModel
{
public:
  CDoubleExp(RooRealVar &m, const Bool_t pass);
  ~CDoubleExp();
  RooExponential *exp1, *exp2;
  RooRealVar *t1, *t2, *frac;
};

class CLinearExp : public CBackgroundModel
{
public:
  CLinearExp(RooRealVar &m, const Bool_t pass);
  ~CLinearExp();
  RooRealVar *a, *t;
};

class CQuadraticExp : public CBackgroundModel
{
public:
  CQuadraticExp(RooRealVar &m, const Bool_t pass);
  ~CQuadraticExp();
  RooRealVar *a1, *a2, *t;
};

//--------------------------------------------------------------------------------------------------
CExponential::CExponential(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  char formula[256];
  
  sprintf(vname,"t%s",name);
  if(pass)
    t = new RooRealVar(vname,vname,-0.04,-1.,1);
  else
    t = new RooRealVar(vname,vname,-0.04,-1.,1);
  sprintf(vname, "expOffset%s", name); offset = new RooRealVar(vname,vname, 0,-200,200);
  sprintf(formula, "m - expOffset%s", name);
  RooArgList *formulaVars = new RooArgList(m, *offset);
  sprintf(vname, "mMinusOffset%s", name); mMinusOffset = new RooFormulaVar(vname, vname, formula, *formulaVars);

      
  sprintf(vname,"background%s",name);
  model = new RooExponential(vname,vname,*mMinusOffset,*t);
}

CExponential::~CExponential()
{
  delete t;
  delete offset;
  delete mMinusOffset;
  t=0;
  offset=0;
}

//--------------------------------------------------------------------------------------------------
CErfcExpo::CErfcExpo(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];

  if(pass) {
    sprintf(vname,"alfa%s",name);  alfa  = new RooRealVar(vname,vname,100,5,200);
    sprintf(vname,"beta%s",name);  beta  = new RooRealVar(vname,vname,0.02,0,0.2);
    sprintf(vname,"gamma%s",name); gamma = new RooRealVar(vname,vname,0.03,0,1);
  } else {
    sprintf(vname,"alfa%s",name);  alfa  = new RooRealVar(vname,vname,100,5,200);
    sprintf(vname,"beta%s",name);  beta  = new RooRealVar(vname,vname,0.02,0,0.2);
    sprintf(vname,"gamma%s",name); gamma = new RooRealVar(vname,vname,0.03,0,1);
  }  
  
  sprintf(vname,"peak%s",name);  
  peak = new RooRealVar(vname,vname,91.1876,85,97); 
  peak->setVal(91.1876);
  peak->setConstant(kTRUE);  
  
  sprintf(vname,"background%s",name);
  model = new RooCMSShape(vname,vname,m,*alfa,*beta,*gamma,*peak);
}

CErfcExpo::~CErfcExpo()
{
  delete alfa;  alfa=0;
  delete beta;  beta=0;
  delete gamma; gamma=0;
  delete peak;  peak=0;
}

//-------------------------------------------------------------------------------------------------
CErfcExpoFixed::CErfcExpoFixed(RooRealVar &m, const Bool_t pass, const int ibin, const std::string fitname, std::string refDir)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  std::vector<std::string> paramNames;
  if(pass) {
    paramNames.push_back("alfaPass");
    paramNames.push_back("betaPass");
    paramNames.push_back("gammaPass");
  } else {
    paramNames.push_back("alfaFail");
    paramNames.push_back("betaFail");
    paramNames.push_back("gammaFail");
  }
  std::vector<double> params = readBkgParams(ibin, fitname, paramNames, refDir);
  sprintf(vname,"alfa%s",name);  alfa  = new RooRealVar(vname,vname,params[0], TMath::Max(5., .1*params[0]), TMath::Min(200., 1.1*params[0]));
  sprintf(vname,"beta%s",name);  beta  = new RooRealVar(vname,vname,params[1]); beta->setConstant(kTRUE);
  sprintf(vname,"gamma%s",name); gamma = new RooRealVar(vname,vname,params[2]); gamma->setConstant(kTRUE);
  
  sprintf(vname,"peak%s",name);  
  peak = new RooRealVar(vname,vname,91.1876,85,97); 
  peak->setVal(91.1876);
  peak->setConstant(kTRUE);  
  
  sprintf(vname,"background%s",name);
  model = new RooCMSShape(vname,vname,m,*alfa,*beta,*gamma,*peak);
}

CErfcExpoFixed::~CErfcExpoFixed()
{
  delete alfa;  alfa=0;
  delete beta;  beta=0;
  delete gamma; gamma=0;
  delete peak;  peak=0;
}

//--------------------------------------------------------------------------------------------------
CDoubleExp::CDoubleExp(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
 
  if(pass) {
    sprintf(vname,"t1%s",name);   t1   = new RooRealVar(vname,vname,-0.20,-1.,0.);
    sprintf(vname,"t2%s",name);   t2   = new RooRealVar(vname,vname,-0.05,-1.,0.);
    sprintf(vname,"frac%s",name); frac = new RooRealVar(vname,vname, 0.50, 0.,1.);
  } else {
    sprintf(vname,"t1%s",name);   t1   = new RooRealVar(vname,vname,-0.20,-1.,0.);
    sprintf(vname,"t2%s",name);   t2   = new RooRealVar(vname,vname,-0.05,-1.,0.);
    sprintf(vname,"frac%s",name); frac = new RooRealVar(vname,vname, 0.50, 0.,1.);
  }
    
  sprintf(vname,"exp1%s",name);
  exp1 = new RooExponential(vname,vname,m,*t1);
  sprintf(vname,"exp2%s",name);
  exp2 = new RooExponential(vname,vname,m,*t2);
  sprintf(vname,"background%s",name);
  model = new RooAddPdf(vname,vname,RooArgList(*exp1,*exp2),RooArgList(*frac));
}

CDoubleExp::~CDoubleExp()
{
  delete exp1; exp1=0;
  delete exp2; exp2=0;
  delete t1;   t1=0;
  delete t2;   t2=0;
  delete frac; frac=0;
}

//--------------------------------------------------------------------------------------------------
CLinearExp::CLinearExp(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char aname[50];
  sprintf(aname,"a%s",name);
  a = new RooRealVar(aname,aname,-0,-10.,10.);
  //a->setConstant(kTRUE);
  
  char tname[50];
  sprintf(tname,"t%s",name);
  t = new RooRealVar(tname,tname,-1e-6,-10.,0.);
  //t->setConstant(kTRUE); 
  
  char formula[200];
  sprintf(formula,"(1+%s*m)*exp(%s*m)",aname,tname);
 
  char vname[50]; sprintf(vname,"background%s",name);
  model = new RooGenericPdf(vname,vname,formula,RooArgList(m,*a,*t));
}

CLinearExp::~CLinearExp()
{
  delete a; a=0;
  delete t; t=0;
}

//--------------------------------------------------------------------------------------------------
CQuadraticExp::CQuadraticExp(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");

  char a1name[50]; 
  sprintf(a1name,"a1%s",name);
  a1 = new RooRealVar(a1name,a1name,0,-10,10.);
  //a1->setConstant(kTRUE);
  
  char a2name[50]; 
  sprintf(a2name,"a2%s",name);
  a2 = new RooRealVar(a2name,a2name,0.0,-10,10);
  //a2->setConstant(kTRUE);
  
  char tname[50];
  sprintf(tname,"t%s",name);
  t = new RooRealVar(tname,tname,-1e-6,-10.,0.); 
  //t->setConstant(kTRUE); 
  
  char formula[200];
  sprintf(formula,"(1+%s*m+%s*m*m)*exp(%s*m)",a1name,a2name,tname);
 
  char vname[50]; sprintf(vname,"background%s",name);
  model = new RooGenericPdf(vname,vname,formula,RooArgList(m,*a1,*a2,*t));
}

CQuadraticExp::~CQuadraticExp()
{
  delete a1; a1=0;
  delete a2; a2=0;
  delete t;  t=0;
}
//--------------------------------------------------------------------------------------------------
std::vector<double> CBackgroundModel::readBkgParams(
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

