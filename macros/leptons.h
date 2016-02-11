#include <TROOT.h>
#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <vector>
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
// LepBaseline = 1UL<<0,
// LepVeto     = 1UL<<1,
// LepFake     = 1UL<<2,
// LepSoft     = 1UL<<3,
// LepLoose    = 1UL<<4,
// LepMedium   = 1UL<<5,
// LepTight    = 1UL<<6,
double selectIsoCut(int type, int pdgId, double eta) {
  bool isEB = TMath::Abs(eta) < 1.479;
  if     (TMath::Abs(pdgId) == 15) {
    if (type==0) return 10000;
    return 4.5;
  }
  if     (TMath::Abs(pdgId) == 13) {
    if (type==0) return 10000;
    return 0.15;
  }
  else if(TMath::Abs(pdgId) == 11) {
    //if     (type == "veto")   return (isEB ? 0.1260 : 0.1440);
    //else if(type == "loose")  return (isEB ? 0.0893 : 0.1210);
    //else if(type == "medium") return (isEB ? 0.0766 : 0.0678);
    //else if(type == "tight")  return (isEB ? 0.0354 : 0.0646);
    if     (type == 0)   return 10000;
    else if(type == 1)   return (isEB ? 0.1260 : 0.1440);
    else if(type == 4)   return (isEB ? 0.0893 : 0.1210);
    else if(type == 5)   return (isEB ? 0.0766 : 0.0678);
    else if(type == 6)   return (isEB ? 0.0354 : 0.0646);
  }
  //printf("Problem with selectIsoCut! type=%d, pdgId=%d, eta=%f\n", type, pdgId,eta);
  //assert(0);
  return 0.0;
}
bool selector(
  int sel_bits,
  double iso,
  double eta,
  double pdgId,
  int id_bit,
  int iso_bit
) {
  if(id_bit < 0) return true;
  if(
    (sel_bits & (0x1 << id_bit)) != 0 &&
    iso < selectIsoCut(iso_bit, pdgId, eta)
  ) return true;
  return false;
}

bool passJetId(Float_t fMVACut[4][4], double mva, double pt, double eta){
                            
  int lPtId = 3;
  if     (pt < 10.)
    lPtId = 0;    
  else if(pt < 20.)
    lPtId = 1;
  else if(pt < 30.)
    lPtId = 2;      
  
  int lEtaId = 3;                                           
  if     (eta < 2.50)
    lEtaId = 0;
  else if(eta < 2.75)
    lEtaId = 1;            
  else if(eta < 3.00)
    lEtaId = 2;
  
  if (mva > fMVACut[lPtId][lEtaId])
    return true;
  
  return false;                         
  
}
void InitializeJetIdCuts(Float_t fMVACut[4][4])
{ 
  float cutValues[4][4] = {
    -0.95, -0.96 ,-0.94, -0.95,
    -0.95, -0.96 ,-0.94, -0.95,
    -0.15, -0.26 ,-0.16, -0.16,
    -0.15, -0.26 ,-0.16, -0.16
  };
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      fMVACut[i][j] = cutValues[i][j];
    }
  }
  
} 

