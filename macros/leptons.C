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
#include <leptons.h>
#include <good_run.h>
using namespace std;

// LepBaseline = 1UL<<0,
// LepVeto     = 1UL<<1,
// LepFake     = 1UL<<2,
// LepSoft     = 1UL<<3,
// LepLoose    = 1UL<<4,
// LepMedium   = 1UL<<5,
// LepTight    = 1UL<<6,

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void make_tnp_skim(
  // a couple extra branches in the tree to know if leptons are true or not
  string list_of_files,
  string output_basename,
  int tag_id            = 6, // tight 
  int tag_iso           = 6, // tight 
  int probe_id          = 0, // baseline
  int probe_iso         = 0, // baseline
  int passing_probe_id  = 6, // tight; use a negative number for a trigger bit
  int passing_probe_iso = 6, // tight
  bool do_electrons = true,
  bool do_muons = false,
  bool do_taus = false,
  bool real_data = false,
  bool verbose = false,
  bool truth_matching = false,
  int electron_trigger=3,
  int muon_trigger=6,
  int tau_trigger=6,
  double tag_pt_min = 30,
  double tag_eta_max = 2.1
) {	
  ifstream ifs(list_of_files.c_str());
  if (!ifs) {
    printf("bad file list\n");
    exit (EXIT_FAILURE);
  }
  string output_dir = "~dhsu/leptonScaleFactors/root/";
  string electron_filename = output_dir +output_basename+"_electronTnP.root";
  string muon_filename     = output_dir +output_basename+"_muonTnP.root";
  string tau_filename      = output_dir +output_basename+"_tauTnP.root";
 
  int min_runNum=99999999, max_runNum=0;
  if(passing_probe_id<0) { electron_trigger=abs(passing_probe_id); muon_trigger=abs(passing_probe_id); }
  //declare output variables
  unsigned int out_runNum, // event ID
  out_lumiSec,
  out_evtNum,
  out_npv, // number of primary vertices
  pass; // whether probe passes requirements
  float        npu=1;                     // mean number of expected pileup
  float        scale1fb=1;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  int          truth_tag, truth_probe;              // tag, probe truth
  TLorentzVector *p4_tag=0, *p4_probe=0;        // tag, probe 4-vector 
  //float        min_probe_gen_deltaR=0;

  TFile *electron_outfile;
  TTree *electron_pair_tree;
  TFile *muon_outfile;
  TTree *muon_pair_tree;
  TFile *tau_outfile;
  TTree *tau_pair_tree;
  
  if(do_electrons) {
    electron_outfile = TFile::Open(electron_filename.c_str(),"RECREATE");
    electron_pair_tree = new TTree("Events", "Electron skim for TnP script");
    electron_pair_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
    electron_pair_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
    electron_pair_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
    electron_pair_tree->Branch("npv",      &out_npv,      "npv/i"      );  
    electron_pair_tree->Branch("pass",     &pass,     "pass/i"     );  
    electron_pair_tree->Branch("npu",      &npu,      "npu/F"      );  
    electron_pair_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
    electron_pair_tree->Branch("mass",     &mass,     "mass/F"     );  
    electron_pair_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
    electron_pair_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
    electron_pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
    electron_pair_tree->Branch("probe", "TLorentzVector", &p4_probe );          
    if(!real_data) {
      electron_pair_tree->Branch("truth_tag",     &truth_tag,     "truth_tag/I"     );  
      electron_pair_tree->Branch("truth_probe",   &truth_probe,   "truth_probe/I"   );  
    }
  }
  if(do_muons) {
    muon_outfile = TFile::Open(muon_filename.c_str(),"RECREATE");
    muon_pair_tree = new TTree("Events", "Muon skim for TnP script");
    muon_pair_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
    muon_pair_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
    muon_pair_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
    muon_pair_tree->Branch("npv",      &out_npv,      "npv/i"      );  
    muon_pair_tree->Branch("pass",     &pass,     "pass/i"     );  
    muon_pair_tree->Branch("npu",      &npu,      "npu/F"      );  
    muon_pair_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
    muon_pair_tree->Branch("mass",     &mass,     "mass/F"     );  
    muon_pair_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
    muon_pair_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
    muon_pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
    muon_pair_tree->Branch("probe", "TLorentzVector", &p4_probe );      
    if(!real_data) {
      muon_pair_tree->Branch("truth_tag",     &truth_tag,     "truth_tag/I"     );  
      muon_pair_tree->Branch("truth_probe",   &truth_probe,   "truth_probe/I"   );  
    }
  }
  if(do_taus) {
    tau_outfile = TFile::Open(tau_filename.c_str(),"RECREATE");
    tau_pair_tree = new TTree("Events", "Muon skim for TnP script");
    tau_pair_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
    tau_pair_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
    tau_pair_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
    tau_pair_tree->Branch("npv",      &out_npv,      "npv/i"      );  
    tau_pair_tree->Branch("pass",     &pass,     "pass/i"     );  
    tau_pair_tree->Branch("npu",      &npu,      "npu/F"      );  
    tau_pair_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
    tau_pair_tree->Branch("mass",     &mass,     "mass/F"     );  
    tau_pair_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
    tau_pair_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
    tau_pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
    tau_pair_tree->Branch("probe", "TLorentzVector", &p4_probe );      
    if(!real_data) {
      tau_pair_tree->Branch("truth_tag",     &truth_tag,     "truth_tag/I"     );  
      tau_pair_tree->Branch("truth_probe",   &truth_probe,   "truth_probe/I"   );  
      //tau_pair_tree->Branch("min_probe_gen_deltaR",   &min_probe_gen_deltaR,   "min_probe_gen_deltaR/F"   );  
    }
  }

  // declare variables to read from nero ntuple
  Int_t           isRealData=0;
  Int_t           runNum=0;
  Int_t           lumiNum=0;
  ULong64_t       eventNum=0;
  Float_t         rho=0;
  TClonesArray    *jetP4=0;
  vector<float>   *jetRawPt=0;
  vector<float>   *jetBdiscr=0;
  vector<float>   *jetBdiscrLegacy=0;
  vector<float>   *jetPuId=0;
  vector<float>   *jetUnc=0;
  vector<float>   *jetQGL=0;
  vector<int>     *jetFlavour=0;
  vector<int>     *jetMatchedPartonPdgId=0;
  vector<int>     *jetMotherPdgId=0;
  vector<int>     *jetGrMotherPdgId=0;
  vector<unsigned int> *jetSelBits=0;
  vector<float>   *jetQ=0;
  vector<float>   *jetQnoPU=0;
  TClonesArray    *jetpuppiP4=0;
  vector<float>   *jetpuppiRawPt=0;
  vector<float>   *jetpuppiBdiscr=0;
  vector<float>   *jetpuppiBdiscrLegacy=0;
  vector<float>   *jetpuppiPuId=0;
  vector<float>   *jetpuppiUnc=0;
  vector<float>   *jetpuppiQGL=0;
  vector<int>     *jetpuppiFlavour=0;
  vector<int>     *jetpuppiMatchedPartonPdgId=0;
  vector<int>     *jetpuppiMotherPdgId=0;
  vector<int>     *jetpuppiGrMotherPdgId=0;
  vector<unsigned int> *jetpuppiSelBits=0;
  vector<float>   *jetpuppiQ=0;
  vector<float>   *jetpuppiQnoPU=0;
  TClonesArray    *fatjetak8P4=0;
  vector<float>   *fatjetak8RawPt=0;
  vector<int>     *fatjetak8Flavour=0;
  vector<float>   *fatjetak8Tau1=0;
  vector<float>   *fatjetak8Tau2=0;
  vector<float>   *fatjetak8Tau3=0;
  vector<float>   *fatjetak8TrimmedMass=0;
  vector<float>   *fatjetak8PrunedMass=0;
  vector<float>   *fatjetak8FilteredMass=0;
  vector<float>   *fatjetak8SoftdropMass=0;
  TClonesArray    *ak8_subjet=0;
  vector<int>     *ak8jet_hasSubjet=0;
  vector<float>   *ak8subjet_btag=0;
  vector<float>   *fatjetak8Hbb=0;
  vector<float>   *fatjetak8topMVA=0;
  TClonesArray    *fatjetca15P4=0;
  vector<float>   *fatjetca15RawPt=0;
  vector<int>     *fatjetca15Flavour=0;
  vector<float>   *fatjetca15Tau1=0;
  vector<float>   *fatjetca15Tau2=0;
  vector<float>   *fatjetca15Tau3=0;
  vector<float>   *fatjetca15TrimmedMass=0;
  vector<float>   *fatjetca15PrunedMass=0;
  vector<float>   *fatjetca15FilteredMass=0;
  vector<float>   *fatjetca15SoftdropMass=0;
  TClonesArray    *ca15_subjet=0;
  vector<int>     *ca15jet_hasSubjet=0;
  vector<float>   *ca15subjet_btag=0;
  vector<float>   *fatjetca15Hbb=0;
  vector<float>   *fatjetca15topMVA=0;
  TClonesArray    *lepP4=0;
  vector<int>     *lepPdgId=0;
  vector<float>   *lepIso=0;
  vector<unsigned int> *lepSelBits=0;
  vector<float>   *lepPfPt=0;
  vector<float>   *lepChIso=0;
  vector<float>   *lepNhIso=0;
  vector<float>   *lepPhoIso=0;
  vector<float>   *lepPuIso=0;
  TClonesArray    *metP4=0;
  Float_t         metSumEtRaw=0;
  vector<float>   *metPtJESUP=0;
  vector<float>   *metPtJESDOWN=0;
  TClonesArray    *metP4_GEN=0;
  TLorentzVector  *metPuppi=0;
  TClonesArray    *metPuppiSyst=0;
  Float_t         metSumEtRawPuppi=0;
  TLorentzVector  *metNoMu=0;
  TLorentzVector  *metNoHF=0;
  Float_t         metSumEtRawNoHF=0;
  TLorentzVector  *pfMet_e3p0=0;
  TLorentzVector  *trackMet=0;
  Float_t         caloMet_Pt=0;
  Float_t         caloMet_Phi=0;
  Float_t         caloMet_SumEt=0;
  TClonesArray    *genP4=0;
  TClonesArray    *genjetP4=0;
  vector<int>     *genPdgId=0;
  Int_t           puTrueInt=0;
  Float_t         mcWeight=0;
  Float_t         pdfQscale=0;
  Float_t         pdfAlphaQED=0;
  Float_t         pdfAlphaQCD=0;
  Float_t         pdfX1=0;
  Float_t         pdfX2=0;
  Int_t           pdfId1=0;
  Int_t           pdfId2=0;
  Float_t         pdfScalePdf=0;
  Float_t         r2f1=0;
  Float_t         r5f1=0;
  Float_t         r1f2=0;
  Float_t         r2f2=0;
  Float_t         r1f5=0;
  Float_t         r5f5=0;
  vector<float>   *pdfRwgt=0;
  TClonesArray    *photonP4=0;
  vector<float>   *photonIso=0;
  vector<float>   *photonSieie=0;
  vector<unsigned int> *photonSelBits=0;
  vector<float>   *photonChIso=0;
  vector<float>   *photonChIsoRC=0;
  vector<float>   *photonNhIso=0;
  vector<float>   *photonNhIsoRC=0;
  vector<float>   *photonPhoIso=0;
  vector<float>   *photonPhoIsoRC=0;
  vector<float>   *photonPuIso=0;
  vector<float>   *photonPuIsoRC=0;
  vector<float>   *photonRawPt=0;
  TClonesArray    *tauP4=0;
  vector<unsigned int> *tauSelBits=0;
  vector<int>     *tauQ=0;
  vector<float>   *tauM=0;
  vector<float>   *tauIso=0;
  vector<float>   *tauChargedIsoPtSum=0;
  vector<float>   *tauNeutralIsoPtSum=0;
  vector<float>   *tauIsoDeltaBetaCorr=0;
  vector<int>     *triggerFired=0;
  vector<float>   *triggerPrescale=0;
  vector<int>     *triggerLeps=0;
  vector<int>     *triggerJets=0;
  vector<int>     *triggerTaus=0;
  vector<int>     *triggerPhotons=0;
  Int_t           npv;
  string input_file_name;
  Long64_t n_events=0;

  TFile *puFile = TFile::Open("~/leptonScaleFactors/puWeights_13TeV_25ns.root","READ");
  TH1D *puWeights = (TH1D*)puFile->Get("puWeights");
  Long64_t totalWeights=0;
  while(getline(ifs,input_file_name)) {

    printf("Opening file \"%s\"\n",input_file_name.c_str());
    TFile *input_file=TFile::Open(input_file_name.c_str(),"READ");
    TTree *events=(TTree*)input_file->Get("nero/events");
    n_events+=events->GetEntriesFast();

    // book branches for nero ntuple
    events->SetBranchAddress("isRealData",&isRealData);
    events->SetBranchAddress("runNum",&runNum);
    events->SetBranchAddress("lumiNum",&lumiNum);
    events->SetBranchAddress("eventNum",&eventNum);
    events->SetBranchAddress("rho",&rho);
    events->SetBranchAddress("jetP4",&jetP4);
    events->SetBranchAddress("jetRawPt",&jetRawPt);
    events->SetBranchAddress("jetBdiscr",&jetBdiscr);
    events->SetBranchAddress("jetBdiscrLegacy",&jetBdiscrLegacy);
    events->SetBranchAddress("jetPuId",&jetPuId);
    events->SetBranchAddress("jetUnc",&jetUnc);
    events->SetBranchAddress("jetQGL",&jetQGL);
    events->SetBranchAddress("jetFlavour",&jetFlavour);
    events->SetBranchAddress("jetMatchedPartonPdgId",&jetMatchedPartonPdgId);
    events->SetBranchAddress("jetMotherPdgId",&jetMotherPdgId);
    events->SetBranchAddress("jetGrMotherPdgId",&jetGrMotherPdgId);
    events->SetBranchAddress("jetQ",&jetQ);
    events->SetBranchAddress("jetQnoPU",&jetQnoPU);
    events->SetBranchAddress("lepP4",&lepP4);
    events->SetBranchAddress("lepPdgId",&lepPdgId);
    events->SetBranchAddress("lepIso",&lepIso);
    events->SetBranchAddress("lepSelBits",&lepSelBits);
    events->SetBranchAddress("lepPfPt",&lepPfPt);
    events->SetBranchAddress("lepChIso",&lepChIso);
    events->SetBranchAddress("lepNhIso",&lepNhIso);
    events->SetBranchAddress("lepPhoIso",&lepPhoIso);
    events->SetBranchAddress("lepPuIso",&lepPuIso);
    events->SetBranchAddress("metP4",&metP4);
    events->SetBranchAddress("metPtJESUP",&metPtJESUP);
    events->SetBranchAddress("metPtJESDOWN",&metPtJESDOWN);
    events->SetBranchAddress("metNoMu",&metNoMu);
    events->SetBranchAddress("pfMet_e3p0",&pfMet_e3p0);
    events->SetBranchAddress("trackMet",&trackMet);
    events->SetBranchAddress("caloMet_Pt",&caloMet_Pt);
    events->SetBranchAddress("caloMet_Phi",&caloMet_Phi);
    events->SetBranchAddress("caloMet_SumEt",&caloMet_SumEt);
    if(!real_data) {
      events->SetBranchAddress("metP4_GEN",&metP4_GEN);
      events->SetBranchAddress("genP4",&genP4);
      events->SetBranchAddress("genjetP4",&genjetP4);
      events->SetBranchAddress("genPdgId",&genPdgId);
      events->SetBranchAddress("puTrueInt",&puTrueInt);
      events->SetBranchAddress("mcWeight",&mcWeight);
    }
    //events->SetBranchAddress("photonP4", &photonP4);
    //events->SetBranchAddress("photonIso", &photonIso);
    //events->SetBranchAddress("photonSieie", &photonSieie);
    //events->SetBranchAddress("photonSelBits", &photonSelBits);
    //events->SetBranchAddress("photonChIso", &photonChIso);
    //events->SetBranchAddress("photonChIsoRC", &photonChIsoRC);
    //events->SetBranchAddress("photonNhIso", &photonNhIso);
    //events->SetBranchAddress("photonNhIsoRC", &photonNhIsoRC);
    //events->SetBranchAddress("photonPhoIso", &photonPhoIso);
    //events->SetBranchAddress("photonPhoIsoRC", &photonPhoIsoRC);
    //events->SetBranchAddress("photonPuIso", &photonPuIso);
    //events->SetBranchAddress("photonPuIsoRC", &photonPuIsoRC);
    //events->SetBranchAddress("photonRawPt", &photonRawPt);
    if(do_taus) {
      events->SetBranchAddress("tauP4", &tauP4);
      events->SetBranchAddress("tauSelBits", &tauSelBits);
      events->SetBranchAddress("tauQ", &tauQ);
      events->SetBranchAddress("tauM", &tauM);
      events->SetBranchAddress("tauIso", &tauIso);
      events->SetBranchAddress("tauChargedIsoPtSum", &tauChargedIsoPtSum);
      events->SetBranchAddress("tauNeutralIsoPtSum", &tauNeutralIsoPtSum);
      events->SetBranchAddress("tauIsoDeltaBetaCorr", &tauIsoDeltaBetaCorr);
    }
    events->SetBranchAddress("triggerFired",&triggerFired);
    events->SetBranchAddress("triggerPrescale",&triggerPrescale);
    events->SetBranchAddress("triggerLeps",&triggerLeps);
    events->SetBranchAddress("triggerJets",&triggerJets);
    events->SetBranchAddress("triggerTaus",&triggerTaus);
    events->SetBranchAddress("triggerPhotons",&triggerPhotons);
    events->SetBranchAddress("npv",&npv);
    
    Long64_t sum_mc_weights=0;
    Long64_t nentries;
    if(!real_data) {
      
      TTree* all_tree=(TTree*)input_file->Get("nero/all");
      all_tree->SetBranchAddress("mcWeight", &mcWeight);
      nentries=all_tree->GetEntries();
      for (Long64_t i=0; i<nentries;i++) {
        all_tree->GetEntry(i);
        sum_mc_weights+=mcWeight / TMath::Abs(mcWeight);
      }
      printf("%lld events in this file\n", sum_mc_weights);
    }

    Long64_t nbytes = 0;
    nentries = events->GetEntries();
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += events->GetEntry(i);
      //printf("runNum %d lumiNum %d\n", runNum, lumiNum);
      if(real_data && !(good_run(runNum, lumiNum) )) { continue; }
      if(real_data && runNum > max_runNum) max_runNum=runNum;
      if(real_data && runNum < min_runNum) min_runNum=runNum;
      //reweighting
      if(!real_data) {
        //if(npv<=30) scale1fb = puWeights->GetBinContent(npv) * mcWeight / sum_mc_weights;
        if(npv<=30) scale1fb = puWeights->GetBinContent(npv) * mcWeight / TMath::Abs(mcWeight);
        else scale1fb=0;
        //printf("event # %lld, scale1fb=%f\b",eventNum,scale1fb);
      }
      unsigned int n_lep = lepPdgId->size();
      unsigned int n_tau;
      if(do_taus) n_tau = tauQ->size();
      //printf("%d leptons in event %lld\n",n_lep,eventNum);

      //###################################################################################################
      // Commented out code to veto events with jets and met for now

      // Get the number of tight leptons
      // if it's 0, there are no tags, so we should nopt waste any more time
      //vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      //for(unsigned int il=0; il<n_lep; il++) {
      //  TLorentzVector *P4=(TLorentzVector*)lepP4->At(il);
      //  if(
      //    selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), 6, 6)
      //  ) {
      //    idTight.push_back(1); idLep.push_back(il); goodIsTight++;
      //  }
      //}
      //if(goodIsTight == 0) continue;

      // Veto events with significant MET
      double minMET  = TMath::Min(((TLorentzVector*)metP4->At(0))->Pt(),(double)trackMet->Pt());
      //if(minMET >= 30) {
      //  if(verbose) printf("rejecting event with MET %f >= 30\n", minMET);
      //  continue;
      //} else {
      //  if(verbose) printf("event passed MET cut with MET %f <= 30\n", minMET); 
      //}
      
      // Veto events with jets
      Float_t fMVACut[4][4];
      InitializeJetIdCuts(fMVACut);

      //vector<int> idJet;
      //bool isBtag = kFALSE;
      //double bDiscrMax = 0.0;
      //double dPhiJetMET = -1.0;
      //for(int nj=0; nj<jetP4->GetEntriesFast(); nj++){
      //  if(((TLorentzVector*)jetP4->At(nj))->Pt() < 10) continue;
      //  bool passId = passJetId(fMVACut, (float)(*jetPuId)[nj], ((TLorentzVector*)(jetP4->At(nj)))->Pt(), TMath::Abs(((TLorentzVector*)(jetP4->At(nj)))->Eta()));
      //  if(((TLorentzVector*)(jetP4->At(nj)))->Pt() < 30) continue;
      //  if(passId == false) continue;        
  
      //  Bool_t isLepton = kFALSE;
      //  // Check delta-R match of jet with the tight leptons
      //  for(unsigned int jl=0; jl < n_lep; jl++) {
      //    if(((TLorentzVector*)(jetP4->At(nj)))->DeltaR( *((TLorentzVector*)(*lepP4)[jl])) < 0.3)
      //      isLepton = kTRUE;
      //  }   
      //  if(isLepton == kTRUE) continue;
  
      //  //if(dPhiJetMET   == -1) dPhiJetMET = TMath::Abs(((TLorentzVector*)(jetP4->At(nj)))->DeltaPhi( *((TLorentzVector*)(*metP4->At(0))) )*180./3.14159265;
  
      //  if(((TLorentzVector*)(jetP4->At(nj)))->Pt() > 15 &&  (float)((*jetBdiscr)[nj]) > bDiscrMax)
      //    bDiscrMax = (float)((*jetBdiscr)[nj]);
  
      //  idJet.push_back(nj);
      //}   
      //if(idJet.size() > 0) {
      //  if(verbose) printf("rejecting event with %lu jets\n", idJet.size());
      //  continue;
      //} else {
      //  if(verbose) printf("event passed Njets cut with %lu jets\n", idJet.size());
      //}
      //###################################################################################################
      
      std::vector<TLorentzVector> p4_ele_tag_, p4_ele_passing_probe_, p4_ele_failing_probe_, p4_mu_tag_, p4_mu_passing_probe_, p4_mu_failing_probe_;
      std::vector<int> q_ele_tag_, q_ele_passing_probe_, q_ele_failing_probe_, q_mu_tag_, q_mu_passing_probe_, q_mu_failing_probe_;
      // record truth info as 1 or 0
      std::vector<int> truth_ele_tag_, truth_ele_passing_probe_, truth_ele_failing_probe_, truth_mu_tag_, truth_mu_passing_probe_, truth_mu_failing_probe_;
      unsigned int il;
      // Loop over the electrons and muons
      if(n_lep != 0 && (do_electrons || do_muons)) { for(il=0; il<n_lep; il++) {
        if(verbose) {
          printf("pdg ID %d:\n",(*lepPdgId)[il] );
          for(int ij = 0; ij<32; ij++) {
            if (((*triggerLeps)[il] & (0x1 << abs(ij))) != 0) printf("passed trigger #%d\n",ij);
            //if (((*triggerFired)[ij] != 0)) printf("passed trigger #%d\n",ij);
          }
        }
        
        // Now find tags and probes
        TLorentzVector *P4=(TLorentzVector*)lepP4->At(il);
        if( P4->Pt() >= 10.) {
          int truth=0;
          int charge=1;
          if((*lepPdgId)[il] > 0)  charge=-1;

          // Loop over gen level info and try to do delta-R match to the lepton
          if(!real_data) {
            unsigned int n_gen = genPdgId->size();
            if(n_gen != 0) { for(unsigned int ig=0; ig<n_gen; ig++) {
              TLorentzVector *gP4 = (TLorentzVector*)genP4->At(ig);
              if(
                ( (*genPdgId)[ig] == (*lepPdgId)[il] ) && 
                ( gP4->DeltaR( *P4) < .1)
              ) truth=1;
            }}
            if(truth==0 && !real_data && truth_matching) continue;
          }
          
          // Record Electron tags and probes
          if(abs( (*lepPdgId)[il]) == 11 && do_electrons) {
            if(
              (passing_probe_iso >= 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_iso < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                ((*triggerLeps)[il] & (0x1 << abs(passing_probe_id))) == 0
              ))
            ) { // make probe selection but fail test selection
              if(verbose) printf("probe failed test ID, relIso = %f, probe iso %f, test iso %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(probe_id, (*lepPdgId)[il], P4->Eta()), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              p4_ele_failing_probe_.push_back((*P4));
              q_ele_failing_probe_.push_back(charge);
              truth_ele_failing_probe_.push_back(truth);
            }
            if(
              (passing_probe_iso >= 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_iso < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                (((*triggerLeps)[il] & (0x1 << abs(passing_probe_id))) != 0)
              ))
            ) { // make probe selection and pass the test selection
              p4_ele_passing_probe_.push_back((*P4));
              if(verbose) printf("passed test ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              q_ele_passing_probe_.push_back(charge);
              truth_ele_passing_probe_.push_back(truth);
            }
            if(
              // additional cuts on tag
              P4->Pt() >= tag_pt_min &&
              TMath::Abs(P4->Eta()) <= tag_eta_max &&
              selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), tag_id, tag_iso) &&
              ( (((*triggerLeps)[il] & (0x1 << electron_trigger)) != 0) || !real_data)
            ) { // pass tag ID, and trigger matching
              p4_ele_tag_.push_back((*P4));
              if(verbose) printf("passed tag ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(tag_id, (*lepPdgId)[il], P4->Eta()));
              q_ele_tag_.push_back(charge);
              truth_ele_tag_.push_back(truth);
            } else {
              if((P4->Pt() < tag_pt_min || TMath::Abs(P4->Eta()) > tag_eta_max)  && verbose) printf("failed tag selection: (pT, eta) = (%f, %f)\n", P4->Pt(), P4->Eta());
              if(! ( (((*triggerLeps)[il] & (0x1 << electron_trigger)) != 0) || !real_data) && verbose) printf("failed tag selection: trigger matching. triggerLeps = %d \n", (*triggerLeps)[il]);
            }
          }
           
          //Record Muon tags and probes
          if(abs( (*lepPdgId)[il] ) == 13 && do_muons) {
            if(
              (passing_probe_iso >= 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_iso < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                ((*triggerLeps)[il] & (0x1 << abs(passing_probe_id))) == 0
              ))
            ) { // make probe selection but fail passing probe selection
              if(verbose) printf("probe failed test ID, relIso = %f, probe iso %f, test iso %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(probe_id, (*lepPdgId)[il], P4->Eta()), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              p4_mu_failing_probe_.push_back((*P4));
              q_mu_failing_probe_.push_back(charge);
              truth_mu_failing_probe_.push_back(truth);
            }
            if(
              (passing_probe_iso >= 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_iso < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                (((*triggerLeps)[il] & (0x1 << abs(passing_probe_id))) != 0)
              ))
            ) { // make probe selection and pass the test selection
              if(verbose) printf("passed test ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              p4_mu_passing_probe_.push_back((*P4));
              q_mu_passing_probe_.push_back(charge);
              truth_mu_passing_probe_.push_back(truth);
            }
          }
          if(abs( (*lepPdgId)[il] ) == 13 && (do_muons)) {
            if(
              P4->Pt() >= tag_pt_min &&
              TMath::Abs(P4->Eta()) <= tag_eta_max &&
              selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), tag_id, tag_iso) &&
              ( (((*triggerLeps)[il] & (0x1 << muon_trigger)) != 0) || !real_data)
            ) { // pass tag muon ID, and trigger matching
              if(verbose) printf("passed tag ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(tag_id, (*lepPdgId)[il], P4->Eta()));
              p4_mu_tag_.push_back((*P4));
              q_mu_tag_.push_back(charge);
              truth_mu_tag_.push_back(truth);
            } else {
              if((P4->Pt() < tag_pt_min || TMath::Abs(P4->Eta()) > tag_eta_max)  && verbose) printf("failed tag selection: (pT, eta) = (%f, %f)\n", P4->Pt(), P4->Eta());
              if(! ( (((*triggerLeps)[il] & (0x1 << muon_trigger)) != 0) || !real_data) && verbose) printf("failed tag selection: trigger matching. triggerLeps = %d \n", (*triggerLeps)[il]);
            }
          }
          
        }
      }}
      
      
      std::vector<TLorentzVector> p4_tau_tag_, p4_tau_passing_probe_, p4_tau_failing_probe_;
      std::vector<int> q_tau_tag_, q_tau_passing_probe_, q_tau_failing_probe_;
      std::vector<int> truth_tau_tag_, truth_tau_passing_probe_, truth_tau_failing_probe_;
      double visibleTauDeltaRCut = 0.1;
      double emVetoDeltaRCut= 0.3;
      // Loop over the taus
      if(n_tau ==0 && verbose && do_taus) printf(" no taus in this event\n");
      
      if(n_tau != 0 && do_taus) { for(unsigned int it=0; it<n_tau; it++) {
        // Now find probes
        TLorentzVector *P4=(TLorentzVector*)tauP4->At(it);
        if( P4->Pt() < 18. || P4->Eta() > 2.3) continue;
        int truth_hadronic=0;
        int charge=(*tauQ)[it];
        bool emVeto=false;
        
        // veto taus that match to a loose electron or muon
        if( ((((*tauSelBits)[it] & (0x1 << 8)) == 0) || ((*tauSelBits)[it] & (0x1 << 10)) == 0) ) {
          if(verbose) printf("\tVetoing event tau with (pT, eta, phi) = (%f, %f, %f) that failed an anti-e or anti-mu discriminator (tauSelBits=%d)\n", P4->Pt(), P4->Eta(), P4->Phi(), (*tauSelBits)[it]);
          emVeto=true;
        }
        il=0;
        while(!emVeto && il<n_lep) {
          TLorentzVector *emP4 = (TLorentzVector*)lepP4->At(il);
          if(
            selector((*lepSelBits)[il], (*lepIso)[il] / emP4->Pt(), emP4->Eta(), abs( (*lepPdgId)[il]), 4,4) &&
            P4->DeltaR(*emP4) < emVetoDeltaRCut
          ) {
            emVeto=true;
            if(verbose) printf(
              "\tVetoing event tau with (pT, eta, phi) = (%f, %f, %f) that passed delta-R cut of %f with a lepton (ID=%d) (%f, %f, %f)\n",
              P4->Pt(),
              P4->Eta(),
              P4->Phi(),
              emVetoDeltaRCut,
              (*lepPdgId)[il],
              emP4->Pt(),
              emP4->Eta(),
              emP4->Phi()
            );
          }
          il++;
        }
        if(emVeto) continue;
        
        if(!real_data) {
          // Loop over gen level info and try to truth match to tau+nu (hadronic)
          if(verbose) printf("\tTrying to truth-match a tau\n");
          unsigned int n_gen = genPdgId->size();
          if(n_gen != 0) { for(unsigned int i_gen=0; i_gen<n_gen; i_gen++) {
            // Try to find a gen tau that could be this particle; if we do, then look for neutrino daughter states
            if ( (*genPdgId)[i_gen] == -15 * charge ) { // Hadronic taus
              TLorentzVector *genTauP4 = (TLorentzVector*)genP4->At(i_gen);
              if(verbose) printf("\t\tFound a gen. tau [%d] with pdg ID %d, (pT, eta, phi) = (%f, %f, %f)\n",i_gen, (*genPdgId)[i_gen], genTauP4->Pt(), genTauP4->Eta(), genTauP4->Phi());
              for(unsigned int i_nu=0; i_nu < n_gen; i_nu++) {
                if( (*genPdgId)[i_nu] != -16 * charge ) continue; // an ID=15 particle has - charge and should decay to tau neutrino ID=16
                TLorentzVector *genNuP4 = (TLorentzVector*)genP4->At(i_nu);
                if(verbose) printf("\t\tFound a gen. nu [%d] with pdg ID %d, (pT, eta, phi) = (%f, %f, %f)\n",i_nu, (*genPdgId)[i_nu], genNuP4->Pt(), genNuP4->Eta(), genNuP4->Phi());
                TLorentzVector visibleGenTauP4 = (*genTauP4) - (*genNuP4); 
                if(verbose) printf("\t\tConstructing visible hadronic gen. tau with (pT, eta, phi) = (%f, %f, %f)\n", visibleGenTauP4.Pt(), visibleGenTauP4.Eta(), visibleGenTauP4.Phi());
                if( visibleGenTauP4.DeltaR( *P4) < visibleTauDeltaRCut) {
                  truth_hadronic=1;
                  if(verbose) printf("\t\tSUCCESS: Event tau with (pT, eta, phi) = (%f, %f, %f) passed delta-R cut of %f with visible hadronic gen. tau!\n", P4->Pt(), P4->Eta(), P4->Phi(), visibleTauDeltaRCut);
                } else {
                  if(verbose) printf("\t\tFAILURE: Event tau with (pT, eta, phi) = (%f, %f, %f) failed delta-R cut of %f with visible hadronic gen. tau!\n", P4->Pt(), P4->Eta(), P4->Phi(), visibleTauDeltaRCut);
                }
              }
            }
          }}
          if(truth_hadronic==0&& !real_data && truth_matching) continue;
        }
        
        // Record Tau probes
        if(
          (passing_probe_id > 0 && (
            selector((*tauSelBits)[it], (*tauIsoDeltaBetaCorr)[it], P4->Eta(), 15, probe_id, probe_iso) &&
            !selector((*tauSelBits)[it], (*tauIsoDeltaBetaCorr)[it], P4->Eta(), 15, passing_probe_id, passing_probe_iso)
          ))
        ) { // make probe selection but fail test selection
          if(verbose) printf("probe failed test ID, tauSelBits = %d, iso = %f, probe iso %f, test iso %f\n", (*tauSelBits)[it], (*tauIsoDeltaBetaCorr)[it], selectIsoCut(probe_id, 15, P4->Eta()), selectIsoCut(passing_probe_id, 15, P4->Eta()));
          p4_tau_failing_probe_.push_back((*P4));
          q_tau_failing_probe_.push_back(charge);
          truth_tau_failing_probe_.push_back(truth_hadronic);
        }
        
        if(!selector((*tauSelBits)[it], (*tauIsoDeltaBetaCorr)[it], P4->Eta(), 15, probe_id, probe_iso) && verbose) printf("failed probe selection, tauSelBits = %d, iso = %f, probe iso %f\n", (*tauSelBits)[it], (*tauIsoDeltaBetaCorr)[it], selectIsoCut(probe_id, 15, P4->Eta()));
        if(
          (passing_probe_id > 0 && (
            selector((*tauSelBits)[it], (*tauIsoDeltaBetaCorr)[it], P4->Eta(), 15, probe_id, probe_iso) &&
            selector((*tauSelBits)[it], (*tauIsoDeltaBetaCorr)[it], P4->Eta(), 15, passing_probe_id, passing_probe_iso)
          ))
          //)) ||
          //(passing_probe_id < 0 && ( // trigger efficiency
          //  selector((*lepSelBits)[it], (*tauIsoDeltaBetaCorr)[it] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[it]), probe_id, probe_iso) &&
          //  (((*triggerLeps)[it] & (0x1 << abs(passing_probe_id))) != 0)
          //))
        ) { // make probe selection and pass the test selection
          p4_tau_passing_probe_.push_back((*P4));
          if(verbose) printf("passed test ID, iso = %f < %f\n",(*tauIsoDeltaBetaCorr)[it], selectIsoCut(passing_probe_id, 15, P4->Eta()));
          q_tau_passing_probe_.push_back(charge);
          truth_tau_passing_probe_.push_back(truth_hadronic);
        }
      }}
      
      // Find Muonic tau tags
      if(n_lep != 0 && (do_taus)) { for(unsigned int il=0; il<n_lep; il++) {
        if(verbose) {
          printf("pdg ID %d:\n",(*lepPdgId)[il] );
          for(int ij = 0; ij<32; ij++) {
            if (((*triggerLeps)[il] & (0x1 << abs(ij))) != 0) printf("passed trigger #%d\n",ij);
          }
        }
        
        TLorentzVector *P4=(TLorentzVector*)lepP4->At(il);
        if(!(
          P4->Pt() >= tag_pt_min &&
          TMath::Abs(P4->Eta()) <= tag_eta_max &&
          (abs( (*lepPdgId)[il] ) != 13)
        )) continue;
        int truth_muonic=0;
        int charge=1;
        if((*lepPdgId)[il] > 0)  charge=-1;

        if(!real_data) {
        // Loop over gen level info and try to truth match to mu+nu+nu (muonic)
          if(verbose) printf("\tTrying to truth-match a tau\n");
          unsigned int n_gen = genPdgId->size();
          if(n_gen != 0) { for(unsigned int i_gen=0; i_gen<n_gen; i_gen++) {
            // Try to find a gen tau that could be this particle; if we do, then look for neutrino daughter states
            if ( (*genPdgId)[i_gen] == -15 * charge ) { // Muonic taus
              TLorentzVector *genTauP4 = (TLorentzVector*)genP4->At(i_gen);
              if(verbose) printf("\t\tFound a gen. tau [%d] with pdg ID %d, (pT, eta, phi) = (%f, %f, %f)\n",i_gen, (*genPdgId)[i_gen], genTauP4->Pt(), genTauP4->Eta(), genTauP4->Phi());
              for(unsigned int i_nu_tau=0; i_nu_tau < n_gen; i_nu_tau++) {
                if( (*genPdgId)[i_nu_tau] != -16 * charge ) continue; // an ID=15 particle has - charge and should decay to tau neutrino ID=16
                TLorentzVector *genNuTauP4 = (TLorentzVector*)genP4->At(i_nu_tau);
                if(verbose) printf("\t\tFound a gen. tau neutrino [%d] with pdg ID %d, (pT, eta, phi) = (%f, %f, %f)\n",i_nu_tau, (*genPdgId)[i_nu_tau], genNuTauP4->Pt(), genNuTauP4->Eta(), genNuTauP4->Phi());
                for(unsigned int i_nu_mu=0; i_nu_mu < n_gen; i_nu_mu++) {
                  if( (*genPdgId)[i_nu_mu] != 14 * charge ) continue; // an ID=15 muonic tau should decay to muon antineutrino ID = -14
                  TLorentzVector *genNuMuP4 = (TLorentzVector*)genP4->At(i_nu_mu);
                  if(verbose) printf("\t\tFound a gen. muon neutrino [%d] with pdg ID %d, (pT, eta, phi) = (%f, %f, %f)\n",i_nu_mu, (*genPdgId)[i_nu_mu], genNuMuP4->Pt(), genNuMuP4->Eta(), genNuMuP4->Phi());
                  TLorentzVector visibleGenTauP4 = (*genTauP4) - (*genNuTauP4) - (*genNuMuP4); 
                  if(verbose) printf("\t\tConstructing visible muonic gen. tau with (pT, eta, phi) = (%f, %f, %f)\n", visibleGenTauP4.Pt(), visibleGenTauP4.Eta(), visibleGenTauP4.Phi());
                  if( visibleGenTauP4.DeltaR( *P4) < visibleTauDeltaRCut) {
                    truth_muonic=1;
                    if(verbose) printf("\t\tSUCCESS: Event muon with (pT, eta, phi) = (%f, %f, %f) passed delta-R cut of %f with visible muonic gen. tau!\n", P4->Pt(), P4->Eta(), P4->Phi(), visibleTauDeltaRCut);
                  } else {
                    if(verbose) printf("\t\tFAILURE: Event muon with (pT, eta, phi) = (%f, %f, %f) failed delta-R cut of %f with visible muonic gen. tau!\n", P4->Pt(), P4->Eta(), P4->Phi(), visibleTauDeltaRCut);
                  }
                }
              }
            }
          }}
          if(truth_muonic==0 && truth_matching) continue;
        }
        
        //Record Muonic Tau tags
        if(
          //P4->Pt() >= 40 &&
          selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), tag_id, tag_iso) &&
          ( (((*triggerLeps)[il] & (0x1 << muon_trigger)) != 0) || !real_data)
        ) { // pass tag muon ID, and trigger matching
          if(verbose) printf("passed tag ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(tag_id, (*lepPdgId)[il], P4->Eta()));
          p4_tau_tag_.push_back((*P4));
          q_tau_tag_.push_back(charge);
          truth_tau_tag_.push_back(truth_muonic);
        } else {
          //if(P4->Pt() < 40 && verbose) printf("failed tag selection: pT = %f < 40\n", P4->Pt());
          if(! ( (((*triggerLeps)[il] & (0x1 << muon_trigger)) != 0) || !real_data) && verbose) printf("failed tag selection: trigger matching. triggerLeps = %d \n", (*triggerLeps)[il]);
        }
        
      }}

      //demote some integers to unsigned
      out_runNum=runNum;
      out_lumiSec=lumiNum;
      out_evtNum=eventNum;
      out_npv=npv;
      npu=puTrueInt;
      
      // associate pairs: electrons
      if(do_electrons) { for(unsigned int iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
        // passing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_passing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_passing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_ele_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a PASSING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            electron_pair_tree->Fill();
          }
        }
        // failing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_failing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_failing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_ele_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( p4_tag->DeltaR(*p4_probe) < .0001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a FAILING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            electron_pair_tree->Fill();
          }
        }

      }}
      
      // associate pairs: muons
      if(do_muons) { for(unsigned int iTag=0; iTag < p4_mu_tag_.size(); iTag++) {
        // passing probes for muons
        for(unsigned int iProbe=0; iProbe < p4_mu_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_mu_tag_[iTag];
          *p4_probe   = p4_mu_passing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_mu_passing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_mu_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( p4_tag->DeltaR(*p4_probe) < .0001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a PASSING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            muon_pair_tree->Fill();
          }
        }
        // failing probes for muons
        for(unsigned int iProbe=0; iProbe < p4_mu_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_mu_tag_[iTag];
          *p4_probe   = p4_mu_failing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_mu_failing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_mu_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( p4_tag->DeltaR(*p4_probe) < .0001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a FAILING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            muon_pair_tree->Fill();
          }
        }

      }}

      // associate pairs: taus (use leptonic tau for tag particle)
      if(do_taus) { for(unsigned int iTag=0; iTag < p4_tau_tag_.size(); iTag++) {
        // passing probes for taus
        for(unsigned int iProbe=0; iProbe < p4_tau_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_tau_tag_[iTag];
          *p4_probe   = p4_tau_passing_probe_[iProbe];
          qtag    = q_tau_tag_[iTag];
          qprobe  = q_tau_passing_probe_[iProbe];
          truth_tag    = truth_tau_tag_[iTag];
          truth_probe  = truth_tau_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( p4_tag->DeltaR(*p4_probe) < .0001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a PASSING tau pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            tau_pair_tree->Fill();
          }
        }
        // failing probes for taus
        for(unsigned int iProbe=0; iProbe < p4_tau_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_tau_tag_[iTag];
          *p4_probe   = p4_tau_failing_probe_[iProbe];
          qtag    = q_tau_tag_[iTag];
          qprobe  = q_tau_failing_probe_[iProbe];
          truth_tag    = truth_tau_tag_[iTag];
          truth_probe  = truth_tau_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( p4_tag->DeltaR(*p4_probe) < .0001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a FAILING tau pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            tau_pair_tree->Fill();
          }
        }

      }}
    }
    input_file->Close();
    totalWeights+=sum_mc_weights;
  }
  //save tnp trees
  if(do_electrons) {
    electron_outfile->cd();
    electron_pair_tree->Write();
    printf("%lld electron pair events\n", electron_pair_tree->GetEntries());
    electron_outfile->Close();
  }
  if(do_muons) {
    muon_outfile->cd();
    muon_pair_tree->Write();
    printf("%lld muon pair events\n", muon_pair_tree->GetEntries());
    muon_outfile->Close();
  }
  if(do_taus) {
    tau_outfile->cd();
    tau_pair_tree->Write();
    printf("%lld tau pair events\n", tau_pair_tree->GetEntries());
    tau_outfile->Close();
  }
  printf("Complete. %lld events processed, total MC events %lld\n", n_events, totalWeights);
  printf("run num [%d,%d] \n", min_runNum, max_runNum); 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*void make_eu_skim(
  // a couple extra branches in the tree to know if leptons are true or not
  string list_of_files,
  string output_basename,
  int tag_id            = 6, // tight 
  int tag_iso           = 6, // tight 
  int probe_id          = 0, // baseline
  int probe_iso         = 0, // baseline
  int passing_probe_id  = 6, // tight; use a negative number for a trigger bit
  int passing_probe_iso = 6, // tight
  bool real_data = false,
  bool verbose = false,
  bool truth_matching = false
) {	
  ifstream ifs(list_of_files.c_str());
  if (!ifs) {
    printf("bad file list\n");
    exit (EXIT_FAILURE);
  }
  string electron_filename="~dhsu/leptonScaleFactors/root/"+output_basename+"_electronTnP.root";
  string muon_filename="~dhsu/leptonScaleFactors/root/"+output_basename+"_muonTnP.root";
 
  int min_runNum=99999999, max_runNum=0;
  int electron_trigger=7, muon_trigger=12;
  if(passing_probe_id<0) { electron_trigger=abs(passing_probe_id); muon_trigger=abs(passing_probe_id); }
  //declare output variables
  unsigned int out_runNum, // event ID
  out_lumiSec,
  out_evtNum,
  out_npv, // number of primary vertices
  pass; // whether probe passes requirements
  float        npu=1;                     // mean number of expected pileup
  float        scale1fb=1;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  int          truth_tag, truth_probe;              // tag, probe truth
  TLorentzVector *p4_tag=0, *p4_probe=0;        // tag, probe 4-vector 
  
  TFile *outfile;
  TTree *pair_tree;
  
  outfile = TFile::Open(electron_filename.c_str(),"RECREATE");
  pair_tree = new TTree("Events", "Electron skim for TnP script");
  pair_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
  pair_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
  pair_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
  pair_tree->Branch("npv",      &out_npv,      "npv/i"      );  
  pair_tree->Branch("pass",     &pass,     "pass/i"     );  
  pair_tree->Branch("npu",      &npu,      "npu/F"      );  
  pair_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
  pair_tree->Branch("mass",     &mass,     "mass/F"     );  
  pair_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
  pair_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
  pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
  pair_tree->Branch("probe", "TLorentzVector", &p4_probe );          
  if(!real_data) {
    pair_tree->Branch("truth_tag",     &truth_tag,     "truth_tag/I"     );  
    pair_tree->Branch("truth_probe",   &truth_probe,   "truth_probe/I"   );  
  }

  // declare variables to read from nero ntuple
  Int_t           isRealData;
  Int_t           runNum;
  Int_t           lumiNum;
  ULong64_t       eventNum;
  Float_t         rho;
  TClonesArray    *jetP4=0;
  vector<float>   *jetRawPt=0;
  vector<float>   *jetBdiscr=0;
  vector<float>   *jetBdiscrLegacy=0;
  vector<float>   *jetPuId=0;
  vector<float>   *jetUnc=0;
  vector<float>   *jetQGL=0;
  vector<int>     *jetFlavour=0;
  vector<int>     *jetMatchedPartonPdgId=0;
  vector<int>     *jetMotherPdgId=0;
  vector<int>     *jetGrMotherPdgId=0;
  vector<bool>    *jetMonojetId=0;
  vector<bool>    *jetMonojetIdLoose=0;
  vector<float>   *jetQ=0;
  vector<float>   *jetQnoPU=0;
  TClonesArray    *lepP4=0;
  vector<int>     *lepPdgId=0;
  vector<float>   *lepIso=0;
  vector<unsigned int> *lepSelBits=0;
  vector<float>   *lepPfPt=0;
  vector<float>   *lepChIso=0;
  vector<float>   *lepNhIso=0;
  vector<float>   *lepPhoIso=0;
  vector<float>   *lepPuIso=0;
  TClonesArray    *metP4=0;
  vector<float>   *metPtJESUP=0;
  vector<float>   *metPtJESDOWN=0;
  TClonesArray    *metP4_GEN=0;
  TLorentzVector  *metNoMu=0;
  TLorentzVector  *pfMet_e3p0=0;
  TLorentzVector  *trackMet=0;
  Float_t         caloMet_Pt;
  Float_t         caloMet_Phi;
  Float_t         caloMet_SumEt;
  TClonesArray    *genP4=0;
  TClonesArray    *genjetP4=0;
  vector<int>     *genPdgId=0;
  Int_t           puTrueInt;
  Float_t         mcWeight;
  Float_t         pdfQscale;
  Float_t         pdfAlphaQED;
  Float_t         pdfAlphaQCD;
  Float_t         pdfX1;
  Float_t         pdfX2;
  Int_t           pdfId1;
  Int_t           pdfId2;
  Float_t         pdfScalePdf;
  TClonesArray    *photonP4=0;
  vector<float>   *photonIso=0;
  vector<float>   *photonSieie=0;
  vector<int>     *photonTightId=0;
  vector<int>     *photonMediumId=0;
  vector<int>     *photonLooseId=0;
  vector<float>   *photonChIso=0;
  vector<float>   *photonChIsoRC=0;
  vector<float>   *photonNhIso=0;
  vector<float>   *photonNhIsoRC=0;
  vector<float>   *photonPhoIso=0;
  vector<float>   *photonPhoIsoRC=0;
  vector<float>   *photonPuIso=0;
  vector<float>   *photonPuIsoRC=0;
  TClonesArray    *tauP4=0;
  vector<float>   *tauId=0;
  vector<int>     *tauQ=0;
  vector<float>   *tauM=0;
  vector<float>   *tauIso=0;
  vector<int>     *triggerFired=0;
  vector<float>   *triggerPrescale=0;
  vector<int>     *triggerLeps=0;
  vector<int>     *triggerJets=0;
  vector<int>     *triggerTaus=0;
  vector<int>     *triggerPhotons=0;
  Int_t           npv;
 
  string input_file_name;
  Long64_t n_events=0;

  TFile *puFile = TFile::Open("~/leptonScaleFactors/puWeights_13TeV_25ns.root","READ");
  TH1D *puWeights = (TH1D*)puFile->Get("puWeights");
  while(getline(ifs,input_file_name)) {

    printf("Opening file \"%s\"\n",input_file_name.c_str());
    TFile *input_file=TFile::Open(input_file_name.c_str(),"READ");
    TTree *events=(TTree*)input_file->Get("nero/events");
    n_events+=events->GetEntriesFast();

    // book branches for nero ntuple
    events->SetBranchAddress("isRealData",&isRealData);
    events->SetBranchAddress("runNum",&runNum);
    events->SetBranchAddress("lumiNum",&lumiNum);
    events->SetBranchAddress("eventNum",&eventNum);
    events->SetBranchAddress("rho",&rho);
    events->SetBranchAddress("jetP4",&jetP4);
    events->SetBranchAddress("jetRawPt",&jetRawPt);
    events->SetBranchAddress("jetBdiscr",&jetBdiscr);
    events->SetBranchAddress("jetBdiscrLegacy",&jetBdiscrLegacy);
    events->SetBranchAddress("jetPuId",&jetPuId);
    events->SetBranchAddress("jetUnc",&jetUnc);
    events->SetBranchAddress("jetQGL",&jetQGL);
    events->SetBranchAddress("jetFlavour",&jetFlavour);
    events->SetBranchAddress("jetMatchedPartonPdgId",&jetMatchedPartonPdgId);
    events->SetBranchAddress("jetMotherPdgId",&jetMotherPdgId);
    events->SetBranchAddress("jetGrMotherPdgId",&jetGrMotherPdgId);
    events->SetBranchAddress("jetQ",&jetQ);
    events->SetBranchAddress("jetQnoPU",&jetQnoPU);
    events->SetBranchAddress("lepP4",&lepP4);
    events->SetBranchAddress("lepPdgId",&lepPdgId);
    events->SetBranchAddress("lepIso",&lepIso);
    events->SetBranchAddress("lepSelBits",&lepSelBits);
    events->SetBranchAddress("lepPfPt",&lepPfPt);
    events->SetBranchAddress("lepChIso",&lepChIso);
    events->SetBranchAddress("lepNhIso",&lepNhIso);
    events->SetBranchAddress("lepPhoIso",&lepPhoIso);
    events->SetBranchAddress("lepPuIso",&lepPuIso);
    events->SetBranchAddress("metP4",&metP4);
    events->SetBranchAddress("metPtJESUP",&metPtJESUP);
    events->SetBranchAddress("metPtJESDOWN",&metPtJESDOWN);
    events->SetBranchAddress("metP4_GEN",&metP4_GEN);
    events->SetBranchAddress("metNoMu",&metNoMu);
    events->SetBranchAddress("pfMet_e3p0",&pfMet_e3p0);
    events->SetBranchAddress("trackMet",&trackMet);
    events->SetBranchAddress("caloMet_Pt",&caloMet_Pt);
    events->SetBranchAddress("caloMet_Phi",&caloMet_Phi);
    events->SetBranchAddress("caloMet_SumEt",&caloMet_SumEt);
    if(!real_data) {
      events->SetBranchAddress("genP4",&genP4);
      events->SetBranchAddress("genjetP4",&genjetP4);
      events->SetBranchAddress("genPdgId",&genPdgId);
      events->SetBranchAddress("puTrueInt",&puTrueInt);
      events->SetBranchAddress("mcWeight",&mcWeight);
    }
    events->SetBranchAddress("triggerFired",&triggerFired);
    events->SetBranchAddress("triggerPrescale",&triggerPrescale);
    events->SetBranchAddress("triggerLeps",&triggerLeps);
    events->SetBranchAddress("triggerJets",&triggerJets);
    events->SetBranchAddress("triggerTaus",&triggerTaus);
    events->SetBranchAddress("triggerPhotons",&triggerPhotons);
    events->SetBranchAddress("npv",&npv);
    
    Float_t sum_mc_weights=1;
    Long64_t nentries = events->GetEntries();
    if(!real_data) {
      sum_mc_weights=0;
      for (Long64_t i=0; i<nentries;i++) {
        events->GetEntry(i);
        sum_mc_weights+=mcWeight;
      }
    }
    Long64_t nbytes = 0;
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += events->GetEntry(i);
      //printf("runNum %d lumiNum %d\n", runNum, lumiNum);
      if(real_data && !(good_run(runNum, lumiNum) )) { continue; }
      if(real_data && runNum > max_runNum) max_runNum=runNum;
      if(real_data && runNum < min_runNum) min_runNum=runNum;
      //reweighting
      if(!real_data) {
        //if(npv<=30) scale1fb = puWeights->GetBinContent(npv) * mcWeight / sum_mc_weights;
        if(npv<=30) scale1fb = puWeights->GetBinContent(npv) * mcWeight / TMath::Abs(mcWeight);
        else scale1fb=0;
        //printf("event # %lld, scale1fb=%f\b",eventNum,scale1fb);
      }
      unsigned int n_lep = lepPdgId->size();
      //printf("%d leptons in event %lld\n",n_lep,eventNum);
      std::vector<TLorentzVector> p4_ele_tag_, p4_ele_passing_probe_, p4_ele_failing_probe_, p4_mu_tag_, p4_mu_passing_probe_, p4_mu_failing_probe_;
      std::vector<int> q_ele_tag_, q_ele_passing_probe_, q_ele_failing_probe_, q_mu_tag_, q_mu_passing_probe_, q_mu_failing_probe_;
      // record truth info as 1 or 0
      std::vector<int> truth_ele_tag_, truth_ele_passing_probe_, truth_ele_failing_probe_, truth_mu_tag_, truth_mu_passing_probe_, truth_mu_failing_probe_;
      // Loop over the leptons
      if(n_lep != 0 && (do_electrons || do_muons)) { for(unsigned int il=0; il<n_lep; il++) {
        if(verbose) {
          printf("pdg ID %d:\n",(*lepPdgId)[il] );
          for(int ij = 0; ij<32; ij++) {
            //if (((*triggerLeps)[il] & (0x1 << abs(ij))) != 0) printf("passed trigger #%d\n",ij);
            if (((*triggerFired)[ij] != 0)) printf("passed trigger #%d\n",ij);
          
          }
        }
        TLorentzVector *P4=(TLorentzVector*)lepP4->At(il);
        if( P4->Pt() >= 10.) {
          int truth=0;
          int charge=-1;
          if((*lepPdgId)[il] > 0)  charge=1;

          // Loop over gen level info and try to do delta-R match to the lepton
          if(!real_data) {
            unsigned int n_gen = genPdgId->size();
            if(n_gen != 0) { for(unsigned int ig=0; ig<n_gen; ig++) {
              TLorentzVector *gP4 = (TLorentzVector*)genP4->At(ig);
              if(
                ( (*genPdgId)[ig] == (*lepPdgId)[il] ) && 
                ( gP4->DeltaR( *P4) < .1)
              ) truth=1;
            }}
            if(truth==0 && !real_data && truth_matching) continue;
          }
          if(abs( (*lepPdgId)[il]) == 11 && do_electrons) {
            if(
              (passing_probe_id > 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_id < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !(((*triggerLeps)[il] & abs(passing_probe_id)) != 0)
              ))
            ) { // make probe selection but fail test selection
              if(verbose) printf("probe failed test ID, relIso = %f, probe iso %f, test iso %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(probe_id, (*lepPdgId)[il], P4->Eta()), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              p4_ele_failing_probe_.push_back((*P4));
              q_ele_failing_probe_.push_back(charge);
              truth_ele_failing_probe_.push_back(truth);
            }
            if(
              (passing_probe_id > 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_id < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                (((*triggerLeps)[il] & abs(passing_probe_id)) != 0)
              ))
            ) { // make probe selection and pass the test selection
              p4_ele_passing_probe_.push_back((*P4));
              if(verbose) printf("passed test ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              q_ele_passing_probe_.push_back(charge);
              truth_ele_passing_probe_.push_back(truth);
            }
            if(
              selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), tag_id, tag_iso) &&
              ( (((*triggerLeps)[il] & electron_trigger) != 0) || !real_data)
            ) { // pass tag ID, and trigger matching
              p4_ele_tag_.push_back((*P4));
              if(verbose) printf("passed tag ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(tag_id, (*lepPdgId)[il], P4->Eta()));
              q_ele_tag_.push_back(charge);
              truth_ele_tag_.push_back(truth);
            }
          }
          if(abs( (*lepPdgId)[il] ) == 13 && do_muons) {
            if(
              (passing_probe_id > 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_id < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !(((*triggerLeps)[il] & abs(passing_probe_id)) != 0)
              ))
            ) { // make probe selection but fail passing probe selection
              if(verbose) printf("probe failed test ID, relIso = %f, probe iso %f, test iso %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(probe_id, (*lepPdgId)[il], P4->Eta()), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              p4_mu_failing_probe_.push_back((*P4));
              q_mu_failing_probe_.push_back(charge);
              truth_mu_failing_probe_.push_back(truth);
            }
            if(
              (passing_probe_id > 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_id < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                (((*triggerLeps)[il] & abs(passing_probe_id)) != 0)
              ))
            ) { // make probe selection and pass the test selection
              if(verbose) printf("passed test ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              p4_mu_passing_probe_.push_back((*P4));
              q_mu_passing_probe_.push_back(charge);
              truth_mu_passing_probe_.push_back(truth);
            }
            if(
              selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), tag_id, tag_iso) &&
              ( (((*triggerLeps)[il] & muon_trigger) != 0) || !real_data)
            ) { // pass tag ele ID, and trigger matching
              if(verbose) printf("passed tag ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(tag_id, (*lepPdgId)[il], P4->Eta()));
              p4_mu_tag_.push_back((*P4));
              q_mu_tag_.push_back(charge);
              truth_mu_tag_.push_back(truth);
            }
          }
          
        }
      }}

      //demote some integers to unsigned
      out_runNum=runNum;
      out_lumiSec=lumiNum;
      out_evtNum=eventNum;
      out_npv=npv;
      npu=puTrueInt;
      
      // associate pairs: electrons
      for(unsigned int iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
        // passing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_passing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_passing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_ele_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a PASSING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            electron_pair_tree->Fill();
          }
        }
        // failing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_failing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_failing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_ele_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a FAILING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            electron_pair_tree->Fill();
          }
        }

      }
      
      // associate pairs: muons
      for(unsigned int iTag=0; iTag < p4_mu_tag_.size(); iTag++) {
        // passing probes for muons
        for(unsigned int iProbe=0; iProbe < p4_mu_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_mu_tag_[iTag];
          *p4_probe   = p4_mu_passing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_mu_passing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_mu_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a PASSING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            muon_pair_tree->Fill();
          }
        }
        // failing probes for muons
        for(unsigned int iProbe=0; iProbe < p4_mu_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_mu_tag_[iTag];
          *p4_probe   = p4_mu_failing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_mu_failing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_mu_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a FAILING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            muon_pair_tree->Fill();
          }
        }

      }
    }
    input_file->Close();
  }
  //save tnp trees
  if(do_electrons) {
    electron_outfile->cd();
    electron_pair_tree->Write();
    printf("%lld electron pair events\n", electron_pair_tree->GetEntries());
    electron_outfile->Close();
  }
  if(do_muons) {
    muon_outfile->cd();
    muon_pair_tree->Write();
    printf("%lld muon pair events\n", muon_pair_tree->GetEntries());
    muon_outfile->Close();
  }
  printf("Complete. %lld events processed\n", n_events);
  printf("run num [%d,%d] \n", min_runNum, max_runNum); 
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*void make_tnp_skim_double_lepton(
  // a couple extra branches in the tree to know if leptons are true or not
  string list_of_files,
  string output_basename,
  int probe_trigger,
  int test_trigger,
  bool real_data = false,
  bool verbose = false,
  bool truth_matching = false
) {	
  ifstream ifs(list_of_files.c_str());
  if (!ifs) {
    printf("bad file list\n");
    exit (EXIT_FAILURE);
  }
  string electron_filename="~dhsu/leptonScaleFactors/root/"+output_basename+"_leptonTnP.root";
  //string electron_filename="~dhsu/leptonScaleFactors/root/"+output_basename+"_electronTnP.root";
  //string muon_filename="~dhsu/leptonScaleFactors/root/"+output_basename+"_muonTnP.root";
 
  int min_runNum=99999999, max_runNum=0;
  int electron_trigger=3, muon_trigger=6; 
  //declare output variables
  unsigned int out_runNum, // event ID
  out_lumiSec,
  out_evtNum,
  out_npv, // number of primary vertices
  pass; // whether probe passes requirements
  float        npu=1;                     // mean number of expected pileup
  float        scale1fb=1;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  int          truth_tag, truth_probe;              // tag, probe truth
  TLorentzVector *p4_tag=0, *p4_probe=0;        // tag, probe 4-vector 
  
  TFile *electron_outfile;
  TTree *electron_pair_tree;
  TFile *muon_outfile;
  TTree *muon_pair_tree;
  
  if(do_electrons) {
    electron_outfile = TFile::Open(electron_filename.c_str(),"RECREATE");
    electron_pair_tree = new TTree("Events", "Electron skim for TnP script");
    electron_pair_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
    electron_pair_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
    electron_pair_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
    electron_pair_tree->Branch("npv",      &out_npv,      "npv/i"      );  
    electron_pair_tree->Branch("pass",     &pass,     "pass/i"     );  
    electron_pair_tree->Branch("npu",      &npu,      "npu/F"      );  
    electron_pair_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
    electron_pair_tree->Branch("mass",     &mass,     "mass/F"     );  
    electron_pair_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
    electron_pair_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
    electron_pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
    electron_pair_tree->Branch("probe", "TLorentzVector", &p4_probe );          
    if(!real_data) {
      electron_pair_tree->Branch("truth_tag",     &truth_tag,     "truth_tag/I"     );  
      electron_pair_tree->Branch("truth_probe",   &truth_probe,   "truth_probe/I"   );  
    }
  }
  if(do_muons) {
    muon_outfile = TFile::Open(muon_filename.c_str(),"RECREATE");
    muon_pair_tree = new TTree("Events", "Muon skim for TnP script");
    muon_pair_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
    muon_pair_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
    muon_pair_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
    muon_pair_tree->Branch("npv",      &out_npv,      "npv/i"      );  
    muon_pair_tree->Branch("pass",     &pass,     "pass/i"     );  
    muon_pair_tree->Branch("npu",      &npu,      "npu/F"      );  
    muon_pair_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
    muon_pair_tree->Branch("mass",     &mass,     "mass/F"     );  
    muon_pair_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
    muon_pair_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
    muon_pair_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
    muon_pair_tree->Branch("probe", "TLorentzVector", &p4_probe );      
    if(!real_data) {
      muon_pair_tree->Branch("truth_tag",     &truth_tag,     "truth_tag/I"     );  
      muon_pair_tree->Branch("truth_probe",   &truth_probe,   "truth_probe/I"   );  
    }
  }

  // declare variables to read from nero ntuple
  Int_t           isRealData;
  Int_t           runNum;
  Int_t           lumiNum;
  ULong64_t       eventNum;
  Float_t         rho;
  TClonesArray    *jetP4=0;
  vector<float>   *jetRawPt=0;
  vector<float>   *jetBdiscr=0;
  vector<float>   *jetBdiscrLegacy=0;
  vector<float>   *jetPuId=0;
  vector<float>   *jetUnc=0;
  vector<float>   *jetQGL=0;
  vector<int>     *jetFlavour=0;
  vector<int>     *jetMatchedPartonPdgId=0;
  vector<int>     *jetMotherPdgId=0;
  vector<int>     *jetGrMotherPdgId=0;
  vector<bool>    *jetMonojetId=0;
  vector<bool>    *jetMonojetIdLoose=0;
  vector<float>   *jetQ=0;
  vector<float>   *jetQnoPU=0;
  TClonesArray    *lepP4=0;
  vector<int>     *lepPdgId=0;
  vector<float>   *lepIso=0;
  vector<unsigned int> *lepSelBits=0;
  vector<float>   *lepPfPt=0;
  vector<float>   *lepChIso=0;
  vector<float>   *lepNhIso=0;
  vector<float>   *lepPhoIso=0;
  vector<float>   *lepPuIso=0;
  TClonesArray    *metP4=0;
  vector<float>   *metPtJESUP=0;
  vector<float>   *metPtJESDOWN=0;
  TClonesArray    *metP4_GEN=0;
  TLorentzVector  *metNoMu=0;
  TLorentzVector  *pfMet_e3p0=0;
  TLorentzVector  *trackMet=0;
  Float_t         caloMet_Pt;
  Float_t         caloMet_Phi;
  Float_t         caloMet_SumEt;
  TClonesArray    *genP4=0;
  TClonesArray    *genjetP4=0;
  vector<int>     *genPdgId=0;
  Int_t           puTrueInt;
  Float_t         mcWeight;
  Float_t         pdfQscale;
  Float_t         pdfAlphaQED;
  Float_t         pdfAlphaQCD;
  Float_t         pdfX1;
  Float_t         pdfX2;
  Int_t           pdfId1;
  Int_t           pdfId2;
  Float_t         pdfScalePdf;
  TClonesArray    *photonP4=0;
  vector<float>   *photonIso=0;
  vector<float>   *photonSieie=0;
  vector<int>     *photonTightId=0;
  vector<int>     *photonMediumId=0;
  vector<int>     *photonLooseId=0;
  vector<float>   *photonChIso=0;
  vector<float>   *photonChIsoRC=0;
  vector<float>   *photonNhIso=0;
  vector<float>   *photonNhIsoRC=0;
  vector<float>   *photonPhoIso=0;
  vector<float>   *photonPhoIsoRC=0;
  vector<float>   *photonPuIso=0;
  vector<float>   *photonPuIsoRC=0;
  TClonesArray    *tauP4=0;
  vector<float>   *tauId=0;
  vector<int>     *tauQ=0;
  vector<float>   *tauM=0;
  vector<float>   *tauIso=0;
  vector<int>     *triggerFired=0;
  vector<float>   *triggerPrescale=0;
  vector<int>     *triggerLeps=0;
  vector<int>     *triggerJets=0;
  vector<int>     *triggerTaus=0;
  vector<int>     *triggerPhotons=0;
  Int_t           npv;
 
  string input_file_name;
  Long64_t n_events=0;

  TFile *puFile = TFile::Open("~/leptonScaleFactors/puWeights_13TeV_25ns.root","READ");
  TH1D *puWeights = (TH1D*)puFile->Get("puWeights");
  while(getline(ifs,input_file_name)) {

    printf("Opening file \"%s\"\n",input_file_name.c_str());
    TFile *input_file=TFile::Open(input_file_name.c_str(),"READ");
    TTree *events=(TTree*)input_file->Get("nero/events");
    n_events+=events->GetEntriesFast();

    // book branches for nero ntuple
    events->SetBranchAddress("isRealData",&isRealData);
    events->SetBranchAddress("runNum",&runNum);
    events->SetBranchAddress("lumiNum",&lumiNum);
    events->SetBranchAddress("eventNum",&eventNum);
    events->SetBranchAddress("rho",&rho);
    events->SetBranchAddress("jetP4",&jetP4);
    events->SetBranchAddress("jetRawPt",&jetRawPt);
    events->SetBranchAddress("jetBdiscr",&jetBdiscr);
    events->SetBranchAddress("jetBdiscrLegacy",&jetBdiscrLegacy);
    events->SetBranchAddress("jetPuId",&jetPuId);
    events->SetBranchAddress("jetUnc",&jetUnc);
    events->SetBranchAddress("jetQGL",&jetQGL);
    events->SetBranchAddress("jetFlavour",&jetFlavour);
    events->SetBranchAddress("jetMatchedPartonPdgId",&jetMatchedPartonPdgId);
    events->SetBranchAddress("jetMotherPdgId",&jetMotherPdgId);
    events->SetBranchAddress("jetGrMotherPdgId",&jetGrMotherPdgId);
    events->SetBranchAddress("jetMonojetId",&jetMonojetId);
    events->SetBranchAddress("jetMonojetIdLoose",&jetMonojetIdLoose);
    events->SetBranchAddress("jetQ",&jetQ);
    events->SetBranchAddress("jetQnoPU",&jetQnoPU);
    events->SetBranchAddress("lepP4",&lepP4);
    events->SetBranchAddress("lepPdgId",&lepPdgId);
    events->SetBranchAddress("lepIso",&lepIso);
    events->SetBranchAddress("lepSelBits",&lepSelBits);
    events->SetBranchAddress("lepPfPt",&lepPfPt);
    events->SetBranchAddress("lepChIso",&lepChIso);
    events->SetBranchAddress("lepNhIso",&lepNhIso);
    events->SetBranchAddress("lepPhoIso",&lepPhoIso);
    events->SetBranchAddress("lepPuIso",&lepPuIso);
    events->SetBranchAddress("metP4",&metP4);
    events->SetBranchAddress("metPtJESUP",&metPtJESUP);
    events->SetBranchAddress("metPtJESDOWN",&metPtJESDOWN);
    events->SetBranchAddress("metP4_GEN",&metP4_GEN);
    events->SetBranchAddress("metNoMu",&metNoMu);
    events->SetBranchAddress("pfMet_e3p0",&pfMet_e3p0);
    events->SetBranchAddress("trackMet",&trackMet);
    events->SetBranchAddress("caloMet_Pt",&caloMet_Pt);
    events->SetBranchAddress("caloMet_Phi",&caloMet_Phi);
    events->SetBranchAddress("caloMet_SumEt",&caloMet_SumEt);
    events->SetBranchAddress("genP4",&genP4);
    events->SetBranchAddress("genjetP4",&genjetP4);
    events->SetBranchAddress("genPdgId",&genPdgId);
    events->SetBranchAddress("puTrueInt",&puTrueInt);
    events->SetBranchAddress("mcWeight",&mcWeight);
    events->SetBranchAddress("pdfQscale",&pdfQscale);
    events->SetBranchAddress("pdfAlphaQED",&pdfAlphaQED);
    events->SetBranchAddress("pdfAlphaQCD",&pdfAlphaQCD);
    events->SetBranchAddress("pdfX1",&pdfX1);
    events->SetBranchAddress("pdfX2",&pdfX2);
    events->SetBranchAddress("pdfId1",&pdfId1);
    events->SetBranchAddress("pdfId2",&pdfId2);
    events->SetBranchAddress("pdfScalePdf",&pdfScalePdf);
    events->SetBranchAddress("photonP4",&photonP4);
    events->SetBranchAddress("photonIso",&photonIso);
    events->SetBranchAddress("photonSieie",&photonSieie);
    events->SetBranchAddress("photonTightId",&photonTightId);
    events->SetBranchAddress("photonMediumId",&photonMediumId);
    events->SetBranchAddress("photonLooseId",&photonLooseId);
    events->SetBranchAddress("photonChIso",&photonChIso);
    events->SetBranchAddress("photonChIsoRC",&photonChIsoRC);
    events->SetBranchAddress("photonNhIso",&photonNhIso);
    events->SetBranchAddress("photonNhIsoRC",&photonNhIsoRC);
    events->SetBranchAddress("photonPhoIso",&photonPhoIso);
    events->SetBranchAddress("photonPhoIsoRC",&photonPhoIsoRC);
    events->SetBranchAddress("photonPuIso",&photonPuIso);
    events->SetBranchAddress("photonPuIsoRC",&photonPuIsoRC);
    events->SetBranchAddress("tauP4",&tauP4);
    events->SetBranchAddress("tauId",&tauId);
    events->SetBranchAddress("tauQ",&tauQ);
    events->SetBranchAddress("tauM",&tauM);
    events->SetBranchAddress("tauIso",&tauIso);
    events->SetBranchAddress("triggerFired",&triggerFired);
    events->SetBranchAddress("triggerPrescale",&triggerPrescale);
    events->SetBranchAddress("triggerLeps",&triggerLeps);
    events->SetBranchAddress("triggerJets",&triggerJets);
    events->SetBranchAddress("triggerTaus",&triggerTaus);
    events->SetBranchAddress("triggerPhotons",&triggerPhotons);
    events->SetBranchAddress("npv",&npv);
    
    Float_t sum_mc_weights=1;
    Long64_t nentries = events->GetEntries();
    if(!real_data) {
      sum_mc_weights=0;
      for (Long64_t i=0; i<nentries;i++) {
        events->GetEntry(i);
        sum_mc_weights+=mcWeight;
      }
    }
    Long64_t nbytes = 0;
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += events->GetEntry(i);
      if(real_data && !(
        
        (runNum == 254231 && ((lumiNum>=1 && lumiNum<=24))) || 
        (runNum == 254232 && ((lumiNum>=1 && lumiNum<=81))) || 
        (runNum == 254790 && ((lumiNum>=90 && lumiNum<=90) || (lumiNum>=93 && lumiNum<=630) || (lumiNum>=633 && lumiNum<=697) || (lumiNum>=701 && lumiNum<=715) || (lumiNum>=719 && lumiNum<=784))) || 
        (runNum == 254852 && ((lumiNum>=47 && lumiNum<=94))) || 
        (runNum == 254879 && ((lumiNum>=52 && lumiNum<=52) || (lumiNum>=54 && lumiNum<=140))) || 
        (runNum == 254906 && ((lumiNum>=1 && lumiNum<=75))) || 
        (runNum == 254907 && ((lumiNum>=1 && lumiNum<=52))) || 
        (runNum == 254914 && ((lumiNum>=32 && lumiNum<=32) || (lumiNum>=34 && lumiNum<=78))) || 
        (runNum == 256630 && ((lumiNum>=5 && lumiNum<=26))) || 
        (runNum == 256673 && ((lumiNum>=55 && lumiNum<=56))) || 
        (runNum == 256674 && ((lumiNum>=1 && lumiNum<=2))) || 
        (runNum == 256675 && ((lumiNum>=1 && lumiNum<=106) || (lumiNum>=111 && lumiNum<=164))) || 
        (runNum == 256676 && ((lumiNum>=1 && lumiNum<=160) || (lumiNum>=162 && lumiNum<=208))) || 
        (runNum == 256677 && ((lumiNum>=1 && lumiNum<=291) || (lumiNum>=293 && lumiNum<=390) || (lumiNum>=392 && lumiNum<=397) || (lumiNum>=400 && lumiNum<=455) || (lumiNum>=457 && lumiNum<=482))) || 
        (runNum == 256801 && ((lumiNum>=73 && lumiNum<=263))) || 
        (runNum == 256842 && ((lumiNum>=131 && lumiNum<=132))) || 
        (runNum == 256843 && ((lumiNum>=1 && lumiNum<=204) || (lumiNum>=207 && lumiNum<=284) || (lumiNum>=286 && lumiNum<=378) || (lumiNum>=380 && lumiNum<=461) || (lumiNum>=463 && lumiNum<=587) || (lumiNum>=598 && lumiNum<=627) || (lumiNum>=630 && lumiNum<=661) || (lumiNum>=1001 && lumiNum<=1034) || (lumiNum>=1036 && lumiNum<=1081) || (lumiNum>=1083 && lumiNum<=1191) || (lumiNum>=1193 && lumiNum<=1193) || (lumiNum>=1195 && lumiNum<=1329) || (lumiNum>=1331 && lumiNum<=1332))) || 
        (runNum == 256866 && ((lumiNum>=34 && lumiNum<=47))) || 
        (runNum == 256867 && ((lumiNum>=1 && lumiNum<=16) || (lumiNum>=19 && lumiNum<=94))) || 
        (runNum == 256868 && ((lumiNum>=5 && lumiNum<=33) || (lumiNum>=35 && lumiNum<=200) || (lumiNum>=202 && lumiNum<=492))) || 
        (runNum == 256869 && ((lumiNum>=1 && lumiNum<=34))) || 
        (runNum == 256926 && ((lumiNum>=35 && lumiNum<=50) || (lumiNum>=53 && lumiNum<=62) || (lumiNum>=64 && lumiNum<=65))) || 
        (runNum == 256941 && ((lumiNum>=1 && lumiNum<=17) || (lumiNum>=19 && lumiNum<=29) || (lumiNum>=103 && lumiNum<=105) || (lumiNum>=107 && lumiNum<=126) || (lumiNum>=129 && lumiNum<=129) || (lumiNum>=131 && lumiNum<=168) || (lumiNum>=170 && lumiNum<=170) || (lumiNum>=175 && lumiNum<=290) || (lumiNum>=293 && lumiNum<=294))) || 
        (runNum == 257461 && ((lumiNum>=44 && lumiNum<=95))) || 
        (runNum == 257531 && ((lumiNum>=5 && lumiNum<=45) || (lumiNum>=50 && lumiNum<=143))) || 
        (runNum == 257599 && ((lumiNum>=42 && lumiNum<=118))) || 
        (runNum == 257613 && ((lumiNum>=14 && lumiNum<=1307))) || 
        (runNum == 257614 && ((lumiNum>=1 && lumiNum<=16))) || 
        (runNum == 257645 && ((lumiNum>=37 && lumiNum<=73) || (lumiNum>=75 && lumiNum<=1096))) || 
        (runNum == 257682 && ((lumiNum>=66 && lumiNum<=366))) || 
        (runNum == 257722 && ((lumiNum>=1 && lumiNum<=19))) || 
        (runNum == 257723 && ((lumiNum>=1 && lumiNum<=1) || (lumiNum>=3 && lumiNum<=108) || (lumiNum>=114 && lumiNum<=148))) || 
        (runNum == 257735 && ((lumiNum>=1 && lumiNum<=15))) || 
        (runNum == 257751 && ((lumiNum>=1 && lumiNum<=463))) || 
        (runNum == 257804 && ((lumiNum>=1 && lumiNum<=17))) || 
        (runNum == 257805 && ((lumiNum>=1 && lumiNum<=249))) || 
        (runNum == 257816 && ((lumiNum>=1 && lumiNum<=385))) || 
        (runNum == 257819 && ((lumiNum>=1 && lumiNum<=248))) || 
        (runNum == 257968 && ((lumiNum>=69 && lumiNum<=326))) || 
        (runNum == 257969 && ((lumiNum>=1 && lumiNum<=634))) || 
        (runNum == 258129 && ((lumiNum>=30 && lumiNum<=124))) || 
        (runNum == 258136 && ((lumiNum>=1 && lumiNum<=60))) || 
        (runNum == 258157 && ((lumiNum>=1 && lumiNum<=56))) || 
        (runNum == 258158 && ((lumiNum>=1 && lumiNum<=1088) || (lumiNum>=1091 && lumiNum<=1786))) || 
        (runNum == 258159 && ((lumiNum>=1 && lumiNum<=501))) || 
        (runNum == 258177 && ((lumiNum>=1 && lumiNum<=342) || (lumiNum>=347 && lumiNum<=724) || (lumiNum>=755 && lumiNum<=1939))) || 
        (runNum == 258211 && ((lumiNum>=43 && lumiNum<=129))) || 
        (runNum == 258213 && ((lumiNum>=1 && lumiNum<=165))) || 
        (runNum == 258214 && ((lumiNum>=1 && lumiNum<=217))) || 
        (runNum == 258215 && ((lumiNum>=1 && lumiNum<=6))) || 
        (runNum == 258287 && ((lumiNum>=45 && lumiNum<=144) || (lumiNum>=148 && lumiNum<=227))) || 
        (runNum == 258403 && ((lumiNum>=1 && lumiNum<=251))) || 
        (runNum == 258425 && ((lumiNum>=3 && lumiNum<=136))) || 
        (runNum == 258426 && ((lumiNum>=1 && lumiNum<=10))) || 
        (runNum == 258427 && ((lumiNum>=1 && lumiNum<=107))) || 
        (runNum == 258428 && ((lumiNum>=1 && lumiNum<=159))) || 
        (runNum == 258432 && ((lumiNum>=1 && lumiNum<=4))) || 
        (runNum == 258434 && ((lumiNum>=1 && lumiNum<=453))) || 
        (runNum == 258440 && ((lumiNum>=1 && lumiNum<=442) || (lumiNum>=444 && lumiNum<=732))) || 
        (runNum == 258444 && ((lumiNum>=1 && lumiNum<=37))) || 
        (runNum == 258445 && ((lumiNum>=1 && lumiNum<=302))) || 
        (runNum == 258446 && ((lumiNum>=1 && lumiNum<=142))) || 
        (runNum == 258448 && ((lumiNum>=2 && lumiNum<=100) || (lumiNum>=102 && lumiNum<=731))) || 
        (runNum == 258655 && ((lumiNum>=60 && lumiNum<=68))) || 
        (runNum == 258656 && ((lumiNum>=1 && lumiNum<=334))) || 
        (runNum == 258694 && ((lumiNum>=23 && lumiNum<=199))) || 
        (runNum == 258702 && ((lumiNum>=52 && lumiNum<=402))) || 
        (runNum == 258703 && ((lumiNum>=1 && lumiNum<=389))) || 
        (runNum == 258705 && ((lumiNum>=1 && lumiNum<=100))) || 
        (runNum == 258706 && ((lumiNum>=1 && lumiNum<=733))) || 
        (runNum == 258712 && ((lumiNum>=1 && lumiNum<=524))) || 
        (runNum == 258713 && ((lumiNum>=1 && lumiNum<=161))) || 
        (runNum == 258714 && ((lumiNum>=1 && lumiNum<=67))) || 
        (runNum == 258741 && ((lumiNum>=22 && lumiNum<=72))) || 
        (runNum == 258742 && ((lumiNum>=2 && lumiNum<=693))) || 
        (runNum == 258745 && ((lumiNum>=1 && lumiNum<=260))) || 
        (runNum == 258749 && ((lumiNum>=1 && lumiNum<=204) || (lumiNum>=220 && lumiNum<=604))) || 
        (runNum == 258750 && ((lumiNum>=1 && lumiNum<=197)))
      )) continue;
      if(real_data && runNum > max_runNum) max_runNum=runNum;
      if(real_data && runNum < min_runNum) min_runNum=runNum;
      //reweighting
      if(!real_data) {
        //if(npv<=30) scale1fb = puWeights->GetBinContent(npv) * mcWeight / sum_mc_weights;
        if(npv<=30) scale1fb = puWeights->GetBinContent(npv) * mcWeight / TMath::Abs(mcWeight);
        else scale1fb=0;
        //printf("event # %lld, scale1fb=%f\b",eventNum,scale1fb);
      }
      unsigned int n_lep = lepPdgId->size();
      std::vector<TLorentzVector> p4_ele_tag_, p4_ele_passing_probe_, p4_ele_failing_probe_, p4_mu_tag_, p4_mu_passing_probe_, p4_mu_failing_probe_;
      std::vector<int> q_ele_tag_, q_ele_passing_probe_, q_ele_failing_probe_, q_mu_tag_, q_mu_passing_probe_, q_mu_failing_probe_;
      // record truth info as 1 or 0
      std::vector<int> truth_ele_tag_, truth_ele_passing_probe_, truth_ele_failing_probe_, truth_mu_tag_, truth_mu_passing_probe_, truth_mu_failing_probe_;
      // Loop over the leptons
      if(n_lep != 0 && (do_electrons || do_muons)) { for(unsigned int il=0; il<n_lep; il++) {
        TLorentzVector *P4=(TLorentzVector*)lepP4->At(il);
        if( P4->Pt() >= 10.) {
          int truth=0;
          int charge=-1;
          if((*lepPdgId)[il] > 0)  charge=1;

          // Loop over gen level info and try to do delta-R match to the lepton
          unsigned int n_gen = genPdgId->size();
          if(n_gen != 0) { for(unsigned int ig=0; ig<n_gen; ig++) {
            TLorentzVector *gP4 = (TLorentzVector*)genP4->At(ig);
            if(
              ( (*genPdgId)[ig] == (*lepPdgId)[il] ) && 
              ( gP4->DeltaR( *P4) < .1)
            ) truth=1;
          }}
          if(truth==0 && !real_data && truth_matching) continue;
          if(abs( (*lepPdgId)[il]) == 11 && do_electrons) {
            if(
              (passing_probe_id > 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_id < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !(((*triggerLeps)[il] & (0x1 << abs(passing_probe_id))) != 0)
              ))
            ) { // make probe selection but fail test selection
              if(verbose) printf("probe failed test ID, relIso = %f, probe iso %f, test iso %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(probe_id, (*lepPdgId)[il], P4->Eta()), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              p4_ele_failing_probe_.push_back((*P4));
              q_ele_failing_probe_.push_back(charge);
              truth_ele_failing_probe_.push_back(truth);
            }
            if(
              (passing_probe_id > 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_id < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                (((*triggerLeps)[il] & (0x1 << abs(passing_probe_id))) != 0)
              ))
            ) { // make probe selection and pass the test selection
              p4_ele_passing_probe_.push_back((*P4));
              if(verbose) printf("passed test ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              q_ele_passing_probe_.push_back(charge);
              truth_ele_passing_probe_.push_back(truth);
            }
            if(
              selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), tag_id, tag_iso) &&
              ((((*triggerLeps)[il] & (0x1 << electron_trigger)) != 0) || !real_data)
            ) { // pass tag ID, and trigger matching
              p4_ele_tag_.push_back((*P4));
              if(verbose) printf("passed tag ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(tag_id, (*lepPdgId)[il], P4->Eta()));
              q_ele_tag_.push_back(charge);
              truth_ele_tag_.push_back(truth);
            }
          }
          if(abs( (*lepPdgId)[il] ) == 13 && do_muons) {
            if(
              (passing_probe_id > 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_id < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                !(((*triggerLeps)[il] & (0x1 << abs(passing_probe_id))) != 0)
              ))
            ) { // make probe selection but fail passing probe selection
              if(verbose) printf("probe failed test ID, relIso = %f, probe iso %f, test iso %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(probe_id, (*lepPdgId)[il], P4->Eta()), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              p4_mu_failing_probe_.push_back((*P4));
              q_mu_failing_probe_.push_back(charge);
              truth_mu_failing_probe_.push_back(truth);
            }
            if(
              (passing_probe_id > 0 && (
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), passing_probe_id, passing_probe_iso)
              )) ||
              (passing_probe_id < 0 && ( // trigger efficiency
                selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), probe_id, probe_iso) &&
                (((*triggerLeps)[il] & (0x1 << abs(passing_probe_id))) != 0)
              ))
            ) { // make probe selection and pass the test selection
              if(verbose) printf("passed test ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(passing_probe_id, (*lepPdgId)[il], P4->Eta()));
              p4_mu_passing_probe_.push_back((*P4));
              q_mu_passing_probe_.push_back(charge);
              truth_mu_passing_probe_.push_back(truth);
            }
            if(
              selector((*lepSelBits)[il], (*lepIso)[il] / P4->Pt(), P4->Eta(), abs( (*lepPdgId)[il]), tag_id, tag_iso) &&
              ( (((*triggerLeps)[il] & (0x1 << muon_trigger)) != 0) || !real_data)
            ) { // pass tag ele ID, and trigger matching
              if(verbose) printf("passed tag ID, relIso = %f < %f\n",(*lepIso)[il] / P4->Pt(), selectIsoCut(tag_id, (*lepPdgId)[il], P4->Eta()));
              p4_mu_tag_.push_back((*P4));
              q_mu_tag_.push_back(charge);
              truth_mu_tag_.push_back(truth);
            }
          }
          
        }
      }}

      //demote some integers to unsigned
      out_runNum=runNum;
      out_lumiSec=lumiNum;
      out_evtNum=eventNum;
      out_npv=npv;
      npu=puTrueInt;
      
      // associate pairs: electrons
      for(unsigned int iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
        // passing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_passing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_passing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_ele_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a PASSING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            electron_pair_tree->Fill();
          }
        }
        // failing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_failing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_failing_probe_[iProbe];
          truth_tag    = truth_ele_tag_[iTag];
          truth_probe  = truth_ele_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a FAILING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            electron_pair_tree->Fill();
          }
        }

      }
      
      // associate pairs: muons
      for(unsigned int iTag=0; iTag < p4_mu_tag_.size(); iTag++) {
        // passing probes for muons
        for(unsigned int iProbe=0; iProbe < p4_mu_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_mu_tag_[iTag];
          *p4_probe   = p4_mu_passing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_mu_passing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_mu_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a PASSING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            muon_pair_tree->Fill();
          }
        }
        // failing probes for muons
        for(unsigned int iProbe=0; iProbe < p4_mu_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_mu_tag_[iTag];
          *p4_probe   = p4_mu_failing_probe_[iProbe];
          qtag    = q_mu_tag_[iTag];
          qprobe  = q_mu_failing_probe_[iProbe];
          truth_tag    = truth_mu_tag_[iTag];
          truth_probe  = truth_mu_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a FAILING muon pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            muon_pair_tree->Fill();
          }
        }

      }
    }
    input_file->Close();
  }
  //save tnp trees
  if(do_electrons) {
    electron_outfile->cd();
    electron_pair_tree->Write();
    printf("%lld electron pair events\n", electron_pair_tree->GetEntries());
    electron_outfile->Close();
  }
  if(do_muons) {
    muon_outfile->cd();
    muon_pair_tree->Write();
    printf("%lld muon pair events\n", muon_pair_tree->GetEntries());
    muon_outfile->Close();
  }
  printf("Complete. %lld events processed\n", n_events);
  printf("run num [%d,%d] \n", min_runNum, max_runNum); 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void correlate_id_iso(
  string list_of_files,
  string output_basename,
  int id_bit           = 6, // tight 
  int iso_bit          = 6, // tight 
  bool do_electrons = true,
  bool do_muons = false
) {	
  Float_t ele_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t ele_eta_bins[] = {0, 1.479, 2.4};
  Int_t n_ele_pt_bins=12;
  Int_t n_ele_eta_bins=2;
  Float_t mu_pt_bins[] = {10,15,20,25,30,40,50,60,70,80,90,100,8000};
  Float_t mu_eta_bins[] = {0, 1.5, 2.1, 2.4};
  Int_t n_mu_pt_bins=12;
  Int_t n_mu_eta_bins=3;

  ifstream ifs(list_of_files.c_str());
  if (!ifs) {
    printf("bad file list\n");
    exit (EXIT_FAILURE);
  }
  string electron_filename=output_basename+"_electronCorrelation.root";
  string muon_filename=output_basename+"_muonCorrelation.root";
 
  TFile *electron_outfile;
  TH2D *ele_corr, *ele_id_sum, *ele_iso_sum, *ele_both_sum;
  TFile *muon_outfile;
  TH2D *mu_corr, *mu_id_sum, *mu_iso_sum, *mu_both_sum;
  if(do_electrons) {
    electron_outfile = TFile::Open(electron_filename.c_str(),"RECREATE");
    ele_corr = new TH2D("ele_corr","Correlation coefficient of ID and isolation VS eta & pT", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
    ele_id_sum   =new TH2D("ele_id_sum","ele_id_sum", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
    ele_iso_sum  =new TH2D("ele_iso_sum","ele_iso_sum", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
    ele_both_sum =new TH2D("ele_both_sum","ele_both_sum", n_ele_eta_bins, ele_eta_bins, n_ele_pt_bins, ele_pt_bins);
  }
  if(do_muons) {
    muon_outfile = TFile::Open(muon_filename.c_str(),"RECREATE");
    mu_corr = new TH2D("mu_corr","Correlation coefficient of ID and isolation VS eta & pT", n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
    mu_id_sum   =new TH2D("mu_id_sum","mu_id_sum", n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
    mu_iso_sum  =new TH2D("mu_iso_sum","mu_iso_sum", n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
    mu_both_sum =new TH2D("mu_both_sum","mu_both_sum", n_mu_eta_bins, mu_eta_bins, n_mu_pt_bins, mu_pt_bins);
  }

  // declare variables to read from nero ntuple
  Int_t           isRealData;
  Int_t           runNum;
  Int_t           lumiNum;
  ULong64_t       eventNum;
  Float_t         rho;
  TClonesArray    *jetP4=0;
  vector<float>   *jetRawPt=0;
  vector<float>   *jetBdiscr=0;
  vector<float>   *jetBdiscrLegacy=0;
  vector<float>   *jetPuId=0;
  vector<float>   *jetUnc=0;
  vector<float>   *jetQGL=0;
  vector<int>     *jetFlavour=0;
  vector<int>     *jetMatchedPartonPdgId=0;
  vector<int>     *jetMotherPdgId=0;
  vector<int>     *jetGrMotherPdgId=0;
  vector<bool>    *jetMonojetId=0;
  vector<bool>    *jetMonojetIdLoose=0;
  vector<float>   *jetQ=0;
  vector<float>   *jetQnoPU=0;
  TClonesArray    *lepP4=0;
  vector<int>     *lepPdgId=0;
  vector<float>   *lepIso=0;
  vector<unsigned int> *lepSelBits=0;
  vector<float>   *lepPfPt=0;
  vector<float>   *lepChIso=0;
  vector<float>   *lepNhIso=0;
  vector<float>   *lepPhoIso=0;
  vector<float>   *lepPuIso=0;
  TClonesArray    *metP4=0;
  vector<float>   *metPtJESUP=0;
  vector<float>   *metPtJESDOWN=0;
  TClonesArray    *metP4_GEN=0;
  TLorentzVector  *metNoMu=0;
  TLorentzVector  *pfMet_e3p0=0;
  TLorentzVector  *trackMet=0;
  Float_t         caloMet_Pt;
  Float_t         caloMet_Phi;
  Float_t         caloMet_SumEt;
  TClonesArray    *genP4=0;
  TClonesArray    *genjetP4=0;
  vector<int>     *genPdgId=0;
  Int_t           puTrueInt;
  Float_t         mcWeight;
  Float_t         pdfQscale;
  Float_t         pdfAlphaQED;
  Float_t         pdfAlphaQCD;
  Float_t         pdfX1;
  Float_t         pdfX2;
  Int_t           pdfId1;
  Int_t           pdfId2;
  Float_t         pdfScalePdf;
  TClonesArray    *photonP4=0;
  vector<float>   *photonIso=0;
  vector<float>   *photonSieie=0;
  vector<int>     *photonTightId=0;
  vector<int>     *photonMediumId=0;
  vector<int>     *photonLooseId=0;
  vector<float>   *photonChIso=0;
  vector<float>   *photonChIsoRC=0;
  vector<float>   *photonNhIso=0;
  vector<float>   *photonNhIsoRC=0;
  vector<float>   *photonPhoIso=0;
  vector<float>   *photonPhoIsoRC=0;
  vector<float>   *photonPuIso=0;
  vector<float>   *photonPuIsoRC=0;
  TClonesArray    *tauP4=0;
  vector<float>   *tauId=0;
  vector<int>     *tauQ=0;
  vector<float>   *tauM=0;
  vector<float>   *tauIso=0;
  vector<int>     *triggerFired=0;
  vector<float>   *triggerPrescale=0;
  vector<int>     *triggerLeps=0;
  vector<int>     *triggerJets=0;
  vector<int>     *triggerTaus=0;
  vector<int>     *triggerPhotons=0;
  Int_t           npv;
 
  string input_file_name;
  Long64_t n_events=0;
  double n_ele=0, n_mu=0;
  while(getline(ifs,input_file_name)) {

    printf("Opening file \"%s\"\n",input_file_name.c_str());
    TFile *input_file=TFile::Open(input_file_name.c_str(),"READ");
    TTree *events=(TTree*)input_file->Get("nero/events");
    n_events+=events->GetEntriesFast();

    // book branches for nero ntuple
    events->SetBranchAddress("isRealData",&isRealData);
    events->SetBranchAddress("runNum",&runNum);
    events->SetBranchAddress("lumiNum",&lumiNum);
    events->SetBranchAddress("eventNum",&eventNum);
    events->SetBranchAddress("rho",&rho);
    events->SetBranchAddress("jetP4",&jetP4);
    events->SetBranchAddress("jetRawPt",&jetRawPt);
    events->SetBranchAddress("jetBdiscr",&jetBdiscr);
    events->SetBranchAddress("jetBdiscrLegacy",&jetBdiscrLegacy);
    events->SetBranchAddress("jetPuId",&jetPuId);
    events->SetBranchAddress("jetUnc",&jetUnc);
    events->SetBranchAddress("jetQGL",&jetQGL);
    events->SetBranchAddress("jetFlavour",&jetFlavour);
    events->SetBranchAddress("jetMatchedPartonPdgId",&jetMatchedPartonPdgId);
    events->SetBranchAddress("jetMotherPdgId",&jetMotherPdgId);
    events->SetBranchAddress("jetGrMotherPdgId",&jetGrMotherPdgId);
    events->SetBranchAddress("jetMonojetId",&jetMonojetId);
    events->SetBranchAddress("jetMonojetIdLoose",&jetMonojetIdLoose);
    events->SetBranchAddress("jetQ",&jetQ);
    events->SetBranchAddress("jetQnoPU",&jetQnoPU);
    events->SetBranchAddress("lepP4",&lepP4);
    events->SetBranchAddress("lepPdgId",&lepPdgId);
    events->SetBranchAddress("lepIso",&lepIso);
    events->SetBranchAddress("lepSelBits",&lepSelBits);
    events->SetBranchAddress("lepPfPt",&lepPfPt);
    events->SetBranchAddress("lepChIso",&lepChIso);
    events->SetBranchAddress("lepNhIso",&lepNhIso);
    events->SetBranchAddress("lepPhoIso",&lepPhoIso);
    events->SetBranchAddress("lepPuIso",&lepPuIso);
    events->SetBranchAddress("metP4",&metP4);
    events->SetBranchAddress("metPtJESUP",&metPtJESUP);
    events->SetBranchAddress("metPtJESDOWN",&metPtJESDOWN);
    events->SetBranchAddress("metP4_GEN",&metP4_GEN);
    events->SetBranchAddress("metNoMu",&metNoMu);
    events->SetBranchAddress("pfMet_e3p0",&pfMet_e3p0);
    events->SetBranchAddress("trackMet",&trackMet);
    events->SetBranchAddress("caloMet_Pt",&caloMet_Pt);
    events->SetBranchAddress("caloMet_Phi",&caloMet_Phi);
    events->SetBranchAddress("caloMet_SumEt",&caloMet_SumEt);
    events->SetBranchAddress("genP4",&genP4);
    events->SetBranchAddress("genjetP4",&genjetP4);
    events->SetBranchAddress("genPdgId",&genPdgId);
    events->SetBranchAddress("puTrueInt",&puTrueInt);
    events->SetBranchAddress("mcWeight",&mcWeight);
    events->SetBranchAddress("pdfQscale",&pdfQscale);
    events->SetBranchAddress("pdfAlphaQED",&pdfAlphaQED);
    events->SetBranchAddress("pdfAlphaQCD",&pdfAlphaQCD);
    events->SetBranchAddress("pdfX1",&pdfX1);
    events->SetBranchAddress("pdfX2",&pdfX2);
    events->SetBranchAddress("pdfId1",&pdfId1);
    events->SetBranchAddress("pdfId2",&pdfId2);
    events->SetBranchAddress("pdfScalePdf",&pdfScalePdf);
    events->SetBranchAddress("photonP4",&photonP4);
    events->SetBranchAddress("photonIso",&photonIso);
    events->SetBranchAddress("photonSieie",&photonSieie);
    events->SetBranchAddress("photonTightId",&photonTightId);
    events->SetBranchAddress("photonMediumId",&photonMediumId);
    events->SetBranchAddress("photonLooseId",&photonLooseId);
    events->SetBranchAddress("photonChIso",&photonChIso);
    events->SetBranchAddress("photonChIsoRC",&photonChIsoRC);
    events->SetBranchAddress("photonNhIso",&photonNhIso);
    events->SetBranchAddress("photonNhIsoRC",&photonNhIsoRC);
    events->SetBranchAddress("photonPhoIso",&photonPhoIso);
    events->SetBranchAddress("photonPhoIsoRC",&photonPhoIsoRC);
    events->SetBranchAddress("photonPuIso",&photonPuIso);
    events->SetBranchAddress("photonPuIsoRC",&photonPuIsoRC);
    events->SetBranchAddress("tauP4",&tauP4);
    events->SetBranchAddress("tauId",&tauId);
    events->SetBranchAddress("tauQ",&tauQ);
    events->SetBranchAddress("tauM",&tauM);
    events->SetBranchAddress("tauIso",&tauIso);
    events->SetBranchAddress("triggerFired",&triggerFired);
    events->SetBranchAddress("triggerPrescale",&triggerPrescale);
    events->SetBranchAddress("triggerLeps",&triggerLeps);
    events->SetBranchAddress("triggerJets",&triggerJets);
    events->SetBranchAddress("triggerTaus",&triggerTaus);
    events->SetBranchAddress("triggerPhotons",&triggerPhotons);
    events->SetBranchAddress("npv",&npv);

    Long64_t nentries = events->GetEntries();
    Long64_t nbytes = 0;
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += events->GetEntry(i);
      
      unsigned int n_lep = lepPdgId->size();
      std::vector<TLorentzVector> p4_ele_tag_, p4_ele_passing_probe_, p4_ele_failing_probe_, p4_mu_tag_, p4_mu_passing_probe_, p4_mu_failing_probe_;
      std::vector<int> q_ele_tag_, q_ele_passing_probe_, q_ele_failing_probe_, q_mu_tag_, q_mu_passing_probe_, q_mu_failing_probe_;
      // Object level loop
      if(n_lep != 0) { for(unsigned int il=0; il<n_lep; il++) {
        TLorentzVector *P4=(TLorentzVector*)lepP4->At(il);
        if( P4->Pt() >= 10.) {
          bool pass_id = ((*lepSelBits)[il] & (0x1 << id_bit) != 0);
          bool pass_iso = (*lepIso)[il] / P4->Pt() < selectIsoCut(iso_bit, abs( (*lepPdgId)[il]), P4->Eta());

          if(abs( (*lepPdgId)[il]) == 11 && do_electrons) {
            if(pass_id) ele_id_sum->Fill(P4->Eta(), P4->Pt());
            if(pass_iso) ele_iso_sum->Fill(P4->Eta(), P4->Pt());
            if(pass_id && pass_iso) ele_both_sum->Fill(P4->Eta(), P4->Pt());
            n_ele++;
          }
          if(abs( (*lepPdgId)[il] ) == 13 && do_muons) {
            if(pass_id) mu_id_sum->Fill(P4->Eta(), P4->Pt());
            if(pass_iso) mu_iso_sum->Fill(P4->Eta(), P4->Pt());
            if(pass_id && pass_iso) mu_both_sum->Fill(P4->Eta(), P4->Pt());
            n_mu++;
          }
        }
      }}
    }
    input_file->Close();
  }
   //save correlation coefficient histos
  if(do_electrons) {
    // compute correlation coefficient for each bin
    for(int i_eta = 1; i_eta <= n_ele_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_ele_pt_bins; i_pt++) {
      int n =  ele_id_sum->GetBin(i_eta, i_pt);
      printf("doing electron bin %d\n",n);
      double v_x, v_y, v_xy; // variance for id, iso, id-iso
      double sum_x = ele_id_sum->GetBinContent(n);
      double sum_y = ele_iso_sum->GetBinContent(n);
      double sum_xy = ele_both_sum->GetBinContent(n);
      v_x = sum_x * (n_ele - sum_x ) / (n_ele*n_ele);
      v_y = sum_y * (n_ele - sum_y ) / (n_ele*n_ele);
      v_xy = sum_xy / n_ele - sum_x * sum_y / (n_ele*n_ele) ;
      ele_corr->SetBinContent(n, v_xy / sqrt(v_x*v_y));
    }}
    electron_outfile->cd();
    ele_corr->Write();
    electron_outfile->Close();
  }
  if(do_muons) {
    // compute correlation coefficient for each bin
    for(int i_eta = 1; i_eta <= n_mu_eta_bins; i_eta++) { for(int i_pt = 1; i_pt <= n_mu_pt_bins; i_pt++) {
      int n =  mu_id_sum->GetBin(i_eta, i_pt);
      printf("doing muon bin %d\n",n);
      double v_x, v_y, v_xy; // variance for id, iso, id-iso
      double sum_x = mu_id_sum->GetBinContent(n);
      double sum_y = mu_iso_sum->GetBinContent(n);
      double sum_xy = mu_both_sum->GetBinContent(n);
      v_x = sum_x * (n_mu - sum_x ) / (n_mu*n_mu);
      v_y = sum_y * (n_mu - sum_y ) / (n_mu*n_mu);
      v_xy = sum_xy / n_mu - sum_x * sum_y / (n_mu*n_mu) ;
      mu_corr->SetBinContent(n, v_xy / sqrt(v_x*v_y));
    }}
    muon_outfile->cd();
    mu_corr->Write();
    muon_outfile->Close();
  }
  printf("Complete. %lld events processed\n", n_events);
  printf("run number [%d, %d]\n", min_runNum, max_runNum);
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void plot_electrons(
  string electron_data_filename = "SingleElectron+Run2015D_electronTnP.root",
  string tt_mc_filename     = "TTJets_electronTnP.root",
  string dy_mc_filename     = "DYJetsToLL_electronTnP.root"
) {
  // Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt
  double tt_scale=15.478*831.76/42394154.;
  double dy_scale=15.478*6025.2/28534018.;
  TCut goodruns = "scale1fb*((runNum == 254231 && ((lumiSec>=1  && lumiSec>=24 ))) || (runNum == 254232 && ((lumiSec>=1  && lumiSec>=81 ))) || (runNum == 254790 && ((lumiSec>=90 && lumiSec>=90) || (lumiSec>=93 && lumiSec<=630) || (lumiSec>=633 && lumiSec<=697) || (lumiSec>=701 && lumiSec<=715) || (lumiSec>=719 && lumiSec<=784))) || (runNum == 254852 && ((lumiSec>=47 && lumiSec>=94))) || (runNum == 254879 && ((lumiSec>=52 && lumiSec>=52) || (lumiSec>=54 && lumiSec<=140))) || (runNum == 254906 && ((lumiSec>=1  && lumiSec>=75))) || (runNum == 254907 && ((lumiSec>=1  && lumiSec>=52))))";

  TFile *electron_data_file = TFile::Open(electron_data_filename.c_str(), "READ");
  TFile *tt_mc_file     = TFile::Open(tt_mc_filename.c_str(), "READ");
  TFile *dy_mc_file     = TFile::Open(dy_mc_filename.c_str(), "READ");
  
  TTree *t_electron_data    = (TTree*) electron_data_file->Get("Events"); 
  TTree *t_tt_mc        = (TTree*) tt_mc_file->Get("Events"); 
  TTree *t_dy_mc        = (TTree*) dy_mc_file->Get("Events"); 

  // plot the mass
  TCanvas *c1= new TCanvas("c1","Invariant dielectron mass");
  TH1F *h_mass_electron_data = new TH1F("h_mass_electron_data","Data",30,60,120);
  TH1F *h_mass_tt_mc = new TH1F("h_mass_tt_mc","TT",30,60,120);
  TH1F *h_mass_dy_mc = new TH1F("h_mass_dy_mc","DY",30,60,120);

  t_electron_data->Draw("mass >> h_mass_electron_data",goodruns);
  t_tt_mc->Draw("mass >> h_mass_tt_mc","scale1fb");
  t_dy_mc->Draw("mass >> h_mass_dy_mc","scale1fb");

  h_mass_tt_mc->Scale(tt_scale);
  h_mass_dy_mc->Scale(dy_scale);
  
  THStack *mass_mc_stack = new THStack("mass_mc_stack","");
  h_mass_tt_mc->SetFillColor(9);
  h_mass_dy_mc->SetFillColor(8);
  h_mass_tt_mc->SetLineColor(1);
  h_mass_dy_mc->SetLineColor(1);
  h_mass_electron_data->SetMarkerColor(1);
  h_mass_electron_data->SetMarkerStyle(8);
  h_mass_electron_data->SetLineColor(1);

  mass_mc_stack->Add(h_mass_tt_mc);
  mass_mc_stack->Add(h_mass_dy_mc);
  c1->SetLogy();
  mass_mc_stack->Draw();
  mass_mc_stack->GetXaxis()->SetTitle("Mass [GeV]");
  mass_mc_stack->GetYaxis()->SetTitle("Events / 2 GeV");
  mass_mc_stack->SetMinimum(.1);
  h_mass_electron_data->Draw("SAME E0 P0");

  TLegend *mass_legend=new TLegend(0.7,0.7,0.85,0.85);
  mass_legend->AddEntry(h_mass_electron_data,"Data","lp");
  mass_legend->AddEntry(h_mass_dy_mc,"DY","f");
  mass_legend->AddEntry(h_mass_tt_mc,"TT","f");
  mass_legend->SetFillColor(0);
  mass_legend->Draw("SAME");
  
  // plot the NPV
  TH1F *h_npv_electron_data = new TH1F("h_npv_electron_data","Data",25,0,50);
  TH1F *h_npv_tt_mc = new TH1F("h_npv_tt_mc","TT",25,0,50);
  TH1F *h_npv_dy_mc = new TH1F("h_npv_dy_mc","DY",25,0,50);
  TCanvas *c2= new TCanvas("c2","Number of primary vertices");
  t_electron_data->Draw("npv >> h_npv_electron_data",goodruns);
  t_tt_mc->Draw("npv >> h_npv_tt_mc","scale1fb");
  t_dy_mc->Draw("npv >> h_npv_dy_mc","scale1fb");
  
  h_npv_tt_mc->Scale(tt_scale);
  h_npv_dy_mc->Scale(dy_scale);
  
  THStack *npv_mc_stack = new THStack("npv_mc_stack","");
  h_npv_tt_mc->SetFillColor(9);
  h_npv_dy_mc->SetFillColor(8);
  h_npv_tt_mc->SetLineColor(1);
  h_npv_dy_mc->SetLineColor(1);
  h_npv_electron_data->SetMarkerColor(1);
  h_npv_electron_data->SetMarkerStyle(8);
  h_npv_electron_data->SetLineColor(1);

  npv_mc_stack->Add(h_npv_tt_mc);
  npv_mc_stack->Add(h_npv_dy_mc);
  //c2->SetLogy();
  npv_mc_stack->Draw();
  npv_mc_stack->GetXaxis()->SetTitle("Primary vertices");
  npv_mc_stack->GetYaxis()->SetTitle("Events / 10 PV");
  h_npv_electron_data->Draw("SAME E0 P0");

  TLegend *npv_legend=new TLegend(0.7,0.7,0.85,0.85);
  npv_legend->AddEntry(h_npv_electron_data,"Data","lp");
  npv_legend->AddEntry(h_npv_dy_mc,"DY","f");
  npv_legend->AddEntry(h_npv_tt_mc,"TT","f");
  npv_legend->SetFillColor(0);
  npv_legend->Draw("SAME");

  // plot ratio of NPV in data to MC
  TH1F *h_npv_mc = new TH1F("h_npv_mc","mc",25,0,50);
  TH1F *h_npv_ratio = new TH1F("h_npv_ratio","ratio",25,0,50);
  *h_npv_mc=(*h_npv_tt_mc)+(*h_npv_dy_mc);
  *h_npv_ratio=(*h_npv_electron_data) / (*h_npv_mc);
  
  TCanvas *c3= new TCanvas("c3","NPV Data / MC");
  
  h_npv_ratio->SetMarkerColor(1);
  h_npv_ratio->SetMarkerStyle(8);
  h_npv_ratio->SetLineColor(1);

  h_npv_ratio->Draw("E0 P0");
  h_npv_ratio->GetXaxis()->SetTitle("Ratio of Primary vertices");
  h_npv_ratio->GetYaxis()->SetTitle("Bins of 10 PV");
  h_npv_ratio->SetMinimum(0);
  h_npv_ratio->SetMaximum(2);
  c3->Update();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void plot_muons(
  string muon_data_filename = "SingleMuon_muonTnP.root",
  string tt_mc_filename     = "TTJets_muonTnP.root",
  string dy_mc_filename     = "DYJetsToLL_muonTnP.root"
) {
  // Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt
  double tt_scale=15.478*831.76/42394154.;
  double dy_scale=15.478*6025.2/28534018.;
  TCut goodruns = "scale1fb*((runNum == 254231 && ((lumiSec>=1  && lumiSec>=24 ))) || (runNum == 254232 && ((lumiSec>=1  && lumiSec>=81 ))) || (runNum == 254790 && ((lumiSec>=90 && lumiSec>=90) || (lumiSec>=93 && lumiSec<=630) || (lumiSec>=633 && lumiSec<=697) || (lumiSec>=701 && lumiSec<=715) || (lumiSec>=719 && lumiSec<=784))) || (runNum == 254852 && ((lumiSec>=47 && lumiSec>=94))) || (runNum == 254879 && ((lumiSec>=52 && lumiSec>=52) || (lumiSec>=54 && lumiSec<=140))) || (runNum == 254906 && ((lumiSec>=1  && lumiSec>=75))) || (runNum == 254907 && ((lumiSec>=1  && lumiSec>=52))))";

  TFile *muon_data_file = TFile::Open(muon_data_filename.c_str(), "READ");
  TFile *tt_mc_file     = TFile::Open(tt_mc_filename.c_str(), "READ");
  TFile *dy_mc_file     = TFile::Open(dy_mc_filename.c_str(), "READ");
  
  TTree *t_muon_data    = (TTree*) muon_data_file->Get("Events"); 
  TTree *t_tt_mc        = (TTree*) tt_mc_file->Get("Events"); 
  TTree *t_dy_mc        = (TTree*) dy_mc_file->Get("Events"); 

  // plot the mass
  TCanvas *c1= new TCanvas("c1","Invariant dimuon mass");
  TH1F *h_mass_muon_data = new TH1F("h_mass_muon_data","Data",30,60,120);
  TH1F *h_mass_tt_mc = new TH1F("h_mass_tt_mc","TT",30,60,120);
  TH1F *h_mass_dy_mc = new TH1F("h_mass_dy_mc","DY",30,60,120);

  t_muon_data->Draw("mass >> h_mass_muon_data",goodruns);
  t_tt_mc->Draw("mass >> h_mass_tt_mc","scale1fb");
  t_dy_mc->Draw("mass >> h_mass_dy_mc","scale1fb");

  h_mass_tt_mc->Scale(tt_scale);
  h_mass_dy_mc->Scale(dy_scale);
  
  THStack *mass_mc_stack = new THStack("mass_mc_stack","");
  h_mass_tt_mc->SetFillColor(9);
  h_mass_dy_mc->SetFillColor(8);
  h_mass_tt_mc->SetLineColor(1);
  h_mass_dy_mc->SetLineColor(1);
  h_mass_muon_data->SetMarkerColor(1);
  h_mass_muon_data->SetMarkerStyle(8);
  h_mass_muon_data->SetLineColor(1);

  mass_mc_stack->Add(h_mass_tt_mc);
  mass_mc_stack->Add(h_mass_dy_mc);
  c1->SetLogy();
  mass_mc_stack->Draw();
  mass_mc_stack->GetXaxis()->SetTitle("Mass [GeV]");
  mass_mc_stack->GetYaxis()->SetTitle("Events / 2 GeV");
  mass_mc_stack->SetMinimum(.1);
  h_mass_muon_data->Draw("SAME E0 P0");

  TLegend *mass_legend=new TLegend(0.7,0.7,0.85,0.85);
  mass_legend->AddEntry(h_mass_muon_data,"Data","lp");
  mass_legend->AddEntry(h_mass_dy_mc,"DY","f");
  mass_legend->AddEntry(h_mass_tt_mc,"TT","f");
  mass_legend->SetFillColor(0);
  mass_legend->Draw("SAME");
  
  // plot the NPV
  TH1F *h_npv_muon_data = new TH1F("h_npv_muon_data","Data",25,0,50);
  TH1F *h_npv_tt_mc = new TH1F("h_npv_tt_mc","TT",25,0,50);
  TH1F *h_npv_dy_mc = new TH1F("h_npv_dy_mc","DY",25,0,50);
  TCanvas *c2= new TCanvas("c2","Number of primary vertices");
  t_muon_data->Draw("npv >> h_npv_muon_data",goodruns);
  t_tt_mc->Draw("npv >> h_npv_tt_mc","scale1fb");
  t_dy_mc->Draw("npv >> h_npv_dy_mc","scale1fb");
  
  h_npv_tt_mc->Scale(tt_scale);
  h_npv_dy_mc->Scale(dy_scale);
  
  THStack *npv_mc_stack = new THStack("npv_mc_stack","");
  h_npv_tt_mc->SetFillColor(9);
  h_npv_dy_mc->SetFillColor(8);
  h_npv_tt_mc->SetLineColor(1);
  h_npv_dy_mc->SetLineColor(1);
  h_npv_muon_data->SetMarkerColor(1);
  h_npv_muon_data->SetMarkerStyle(8);
  h_npv_muon_data->SetLineColor(1);

  npv_mc_stack->Add(h_npv_tt_mc);
  npv_mc_stack->Add(h_npv_dy_mc);
  //c2->SetLogy();
  npv_mc_stack->Draw();
  npv_mc_stack->GetXaxis()->SetTitle("Primary vertices");
  npv_mc_stack->GetYaxis()->SetTitle("Events / 10 PV");
  h_npv_muon_data->Draw("SAME E0 P0");

  TLegend *npv_legend=new TLegend(0.7,0.7,0.85,0.85);
  npv_legend->AddEntry(h_npv_muon_data,"Data","lp");
  npv_legend->AddEntry(h_npv_dy_mc,"DY","f");
  npv_legend->AddEntry(h_npv_tt_mc,"TT","f");
  npv_legend->SetFillColor(0);
  npv_legend->Draw("SAME");

  // plot ratio of NPV in data to MC
  TH1F *h_npv_mc = new TH1F("h_npv_mc","mc",25,0,50);
  TH1F *h_npv_ratio = new TH1F("h_npv_ratio","ratio",25,0,50);
  *h_npv_mc=(*h_npv_tt_mc)+(*h_npv_dy_mc);
  *h_npv_ratio=(*h_npv_muon_data) / (*h_npv_mc);
  
  TCanvas *c3= new TCanvas("c3","NPV Data / MC");
  
  h_npv_ratio->SetMarkerColor(1);
  h_npv_ratio->SetMarkerStyle(8);
  h_npv_ratio->SetLineColor(1);

  //c2->SetLogy();
  h_npv_ratio->Draw("E0 P0");
  h_npv_ratio->GetXaxis()->SetTitle("Ratio of Primary vertices");
  h_npv_ratio->GetYaxis()->SetTitle("Bins of 10 PV");
  h_npv_ratio->SetMinimum(0);
  h_npv_ratio->SetMaximum(2);
  c3->Update();

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void reweighted_electron_skim(
  string output_basename,
  int tag_id           = 5, // tight 
  int probe_id         = 3, // loose
  int passing_probe_id = 5,  // tight
  bool verbose=false
) {
  // First read existing skims to compute ratio of npv
  double tt_scale=15.478*831.76/42394154.;
  double dy_scale=15.478*6025.2/28534018.;
  TCut goodruns="";
  
  // loop over MC and data NTUPLE files and build up NPV distributions for each
  // then take the quotient as another histogram which we use to make the reweighted TNP Skim

  ifstream tt_ifs("TTJets_files.txt");
  ifstream dy_ifs("DYJetsToLL_files.txt");
  ifstream data_ifs("SingleElectron_files.txt");

  if (!tt_ifs || !dy_ifs || !data_ifs) {
    if(!tt_ifs) printf("bad tt file list\n");
    if(!dy_ifs) printf("bad dy file list\n");
    if(!data_ifs) printf("bad data file list\n");
    exit (EXIT_FAILURE);
  }
  string input_file_name;
  
  TH1F *h_npv_tt_mc_sum = new TH1F("h_npv_tt_mc","TT",50,0,50);
  TH1F *h_npv_dy_mc_sum = new TH1F("h_npv_dy_mc","DY",50,0,50);
  TH1F *h_npv_electron_data_sum = new TH1F("h_npv_electron_data","Data",50,0,50);
  printf(" Reading the file lists to calculate NPV ratio.\n");
  while(getline(tt_ifs,input_file_name)) {
    if(verbose) printf("Opening file \"%s\"\n",input_file_name.c_str());
    TFile *input_file=TFile::Open(input_file_name.c_str(),"READ");
    TTree *events=(TTree*)input_file->Get("nero/events");
    TH1F *h_npv_tt_mc = new TH1F("h_npv_tt_mc","TT",50,0,50);
    events->Draw("npv >> h_npv_tt_mc");
    h_npv_tt_mc_sum->Add(h_npv_tt_mc);
    input_file->Close();
  }
  while(getline(dy_ifs,input_file_name)) {
    if(verbose) printf("Opening file \"%s\"\n",input_file_name.c_str());
    TFile *input_file=TFile::Open(input_file_name.c_str(),"READ");
    TTree *events=(TTree*)input_file->Get("nero/events");
    TH1F *h_npv_dy_mc = new TH1F("h_npv_dy_mc","DY",50,0,50);
    events->Draw("npv >> h_npv_dy_mc");
    h_npv_dy_mc_sum->Add(h_npv_dy_mc);
    input_file->Close();
  }
  while(getline(data_ifs,input_file_name)) {
    if(verbose) printf("Opening file \"%s\"\n",input_file_name.c_str());
    TFile *input_file=TFile::Open(input_file_name.c_str(),"READ");
    TTree *events=(TTree*)input_file->Get("nero/events");
    TH1F *h_npv_data_mc = new TH1F("h_npv_electron_data","data",50,0,50);
    events->Draw("npv >> h_npv_electron_data",goodruns);
    h_npv_electron_data_sum->Add(h_npv_data_mc);
    input_file->Close();
  }
  tt_ifs.clear();    tt_ifs.seekg(0, std::ios::beg); 
  dy_ifs.clear();    dy_ifs.seekg(0, std::ios::beg); 
  data_ifs.clear();  data_ifs.seekg(0, std::ios::beg); 

  // get the NPV
  
  h_npv_tt_mc_sum->Scale(tt_scale);
  h_npv_dy_mc_sum->Scale(dy_scale);
  
  // get ratio of NPV in data to MC
  TH1F *h_npv_mc = new TH1F("h_npv_mc","mc",50,0,50);
  TH1F *h_npv_ratio = new TH1F("h_npv_ratio","ratio",50,0,50);
  *h_npv_mc=(*h_npv_tt_mc_sum)+(*h_npv_dy_mc_sum);
  *h_npv_ratio=(*h_npv_electron_data_sum) / (*h_npv_mc);
  if(verbose) { for(int j = 1; j <= 50; j++) {
    printf("\tnpv=%d: data=%f, mc=%f, ratio=%f \n", j, h_npv_electron_data_sum->GetBinContent(j), h_npv_mc->GetBinContent(j), h_npv_ratio->GetBinContent(j));
  }}
  //declare output variables
  unsigned int out_runNum, // event ID
  out_lumiSec,
  out_evtNum,
  out_npv, // number of primary vertices
  pass; // whether probe passes requirements
  float        npu=0;                     // mean number of expected pileup
  float        scale1fb=1;                  // event weight per 1/fb
  float        mass;                      // tag-probe mass
  int          qtag, qprobe;              // tag, probe charge
  TLorentzVector *p4_tag=0, *p4_probe=0;        // tag, probe 4-vector 

  //create reweighted files and trees
  string reweighted_tt_filename = "TTJets_"+output_basename+"_electronTnP_reweightedPU.root";
  TFile *reweighted_tt_file     = TFile::Open(reweighted_tt_filename.c_str(),"RECREATE");
  TTree *reweighted_tt_tree = new TTree("Events", "Reweighted TT skim for TnP script");
  string reweighted_dy_filename = "DYJetsToLL_"+output_basename+"_electronTnP_reweightedPU.root";
  TFile *reweighted_dy_file     = TFile::Open(reweighted_dy_filename.c_str(),"RECREATE");
  TTree *reweighted_dy_tree = new TTree("Events", "Reweighted DY skim for TnP script");
  
  //book branches for reweighted trees
  reweighted_tt_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
  reweighted_tt_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
  reweighted_tt_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
  reweighted_tt_tree->Branch("npv",      &out_npv,      "npv/i"      );  
  reweighted_tt_tree->Branch("pass",     &pass,     "pass/i"     );  
  reweighted_tt_tree->Branch("npu",      &npu,      "npu/F"      );  
  reweighted_tt_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
  reweighted_tt_tree->Branch("mass",     &mass,     "mass/F"     );  
  reweighted_tt_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
  reweighted_tt_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
  reweighted_tt_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
  reweighted_tt_tree->Branch("probe", "TLorentzVector", &p4_probe );      
  reweighted_dy_tree->Branch("runNum",   &out_runNum,   "runNum/i"   );  
  reweighted_dy_tree->Branch("lumiSec",  &out_lumiSec,  "lumiSec/i"  );  
  reweighted_dy_tree->Branch("evtNum",   &out_evtNum,   "evtNum/i"   );  
  reweighted_dy_tree->Branch("npv",      &out_npv,      "npv/i"      );  
  reweighted_dy_tree->Branch("pass",     &pass,     "pass/i"     );  
  reweighted_dy_tree->Branch("npu",      &npu,      "npu/F"      );  
  reweighted_dy_tree->Branch("scale1fb", &scale1fb, "scale1fb/F" );
  reweighted_dy_tree->Branch("mass",     &mass,     "mass/F"     );  
  reweighted_dy_tree->Branch("qtag",     &qtag,     "qtag/I"     );  
  reweighted_dy_tree->Branch("qprobe",   &qprobe,   "qprobe/I"   );  
  reweighted_dy_tree->Branch("tag",   "TLorentzVector", &p4_tag   );  
  reweighted_dy_tree->Branch("probe", "TLorentzVector", &p4_probe );      
  
  // declare input variables to read from tree
  Int_t           isRealData;
  Int_t           runNum;
  Int_t           lumiNum;
  ULong64_t       eventNum;
  Float_t         rho;
  TClonesArray    *jetP4=0;
  vector<float>   *jetRawPt=0;
  vector<float>   *jetBdiscr=0;
  vector<float>   *jetBdiscrLegacy=0;
  vector<float>   *jetPuId=0;
  vector<float>   *jetUnc=0;
  vector<float>   *jetQGL=0;
  vector<int>     *jetFlavour=0;
  vector<int>     *jetMatchedPartonPdgId=0;
  vector<int>     *jetMotherPdgId=0;
  vector<int>     *jetGrMotherPdgId=0;
  vector<bool>    *jetMonojetId=0;
  vector<bool>    *jetMonojetIdLoose=0;
  vector<float>   *jetQ=0;
  vector<float>   *jetQnoPU=0;
  TClonesArray    *lepP4=0;
  vector<int>     *lepPdgId=0;
  vector<float>   *lepIso=0;
  vector<unsigned int> *lepSelBits=0;
  vector<float>   *lepPfPt=0;
  vector<float>   *lepChIso=0;
  vector<float>   *lepNhIso=0;
  vector<float>   *lepPhoIso=0;
  vector<float>   *lepPuIso=0;
  TClonesArray    *metP4=0;
  vector<float>   *metPtJESUP=0;
  vector<float>   *metPtJESDOWN=0;
  TClonesArray    *metP4_GEN=0;
  TLorentzVector  *metNoMu=0;
  TLorentzVector  *pfMet_e3p0=0;
  TLorentzVector  *trackMet=0;
  Float_t         caloMet_Pt;
  Float_t         caloMet_Phi;
  Float_t         caloMet_SumEt;
  TClonesArray    *genP4=0;
  TClonesArray    *genjetP4=0;
  vector<int>     *genPdgId=0;
  Int_t           puTrueInt;
  Float_t         mcWeight;
  Float_t         pdfQscale;
  Float_t         pdfAlphaQED;
  Float_t         pdfAlphaQCD;
  Float_t         pdfX1;
  Float_t         pdfX2;
  Int_t           pdfId1;
  Int_t           pdfId2;
  Float_t         pdfScalePdf;
  TClonesArray    *photonP4=0;
  vector<float>   *photonIso=0;
  vector<float>   *photonSieie=0;
  vector<int>     *photonTightId=0;
  vector<int>     *photonMediumId=0;
  vector<int>     *photonLooseId=0;
  vector<float>   *photonChIso=0;
  vector<float>   *photonChIsoRC=0;
  vector<float>   *photonNhIso=0;
  vector<float>   *photonNhIsoRC=0;
  vector<float>   *photonPhoIso=0;
  vector<float>   *photonPhoIsoRC=0;
  vector<float>   *photonPuIso=0;
  vector<float>   *photonPuIsoRC=0;
  TClonesArray    *tauP4=0;
  vector<float>   *tauId=0;
  vector<int>     *tauQ=0;
  vector<float>   *tauM=0;
  vector<float>   *tauIso=0;
  vector<int>     *triggerFired=0;
  vector<float>   *triggerPrescale=0;
  vector<int>     *triggerLeps=0;
  vector<int>     *triggerJets=0;
  vector<int>     *triggerTaus=0;
  vector<int>     *triggerPhotons=0;
  Int_t           npv;
    
  // make TT skim
  printf("making TT skim\n");
  while(getline(tt_ifs,input_file_name)) {
    if(verbose) printf("Opening file \"%s\"\n",input_file_name.c_str());
    TFile *input_file=TFile::Open(input_file_name.c_str(),"READ");
    TTree *events=(TTree*)input_file->Get("nero/events");
    // Set branch addresses.
    events->SetBranchAddress("isRealData",&isRealData);
    events->SetBranchAddress("runNum",&runNum);
    events->SetBranchAddress("lumiNum",&lumiNum);
    events->SetBranchAddress("eventNum",&eventNum);
    events->SetBranchAddress("rho",&rho);
    events->SetBranchAddress("jetP4",&jetP4);
    events->SetBranchAddress("jetRawPt",&jetRawPt);
    events->SetBranchAddress("jetBdiscr",&jetBdiscr);
    events->SetBranchAddress("jetBdiscrLegacy",&jetBdiscrLegacy);
    events->SetBranchAddress("jetPuId",&jetPuId);
    events->SetBranchAddress("jetUnc",&jetUnc);
    events->SetBranchAddress("jetQGL",&jetQGL);
    events->SetBranchAddress("jetFlavour",&jetFlavour);
    events->SetBranchAddress("jetMatchedPartonPdgId",&jetMatchedPartonPdgId);
    events->SetBranchAddress("jetMotherPdgId",&jetMotherPdgId);
    events->SetBranchAddress("jetGrMotherPdgId",&jetGrMotherPdgId);
    events->SetBranchAddress("jetMonojetId",&jetMonojetId);
    events->SetBranchAddress("jetMonojetIdLoose",&jetMonojetIdLoose);
    events->SetBranchAddress("jetQ",&jetQ);
    events->SetBranchAddress("jetQnoPU",&jetQnoPU);
    events->SetBranchAddress("lepP4",&lepP4);
    events->SetBranchAddress("lepPdgId",&lepPdgId);
    events->SetBranchAddress("lepIso",&lepIso);
    events->SetBranchAddress("lepSelBits",&lepSelBits);
    events->SetBranchAddress("lepPfPt",&lepPfPt);
    events->SetBranchAddress("lepChIso",&lepChIso);
    events->SetBranchAddress("lepNhIso",&lepNhIso);
    events->SetBranchAddress("lepPhoIso",&lepPhoIso);
    events->SetBranchAddress("lepPuIso",&lepPuIso);
    events->SetBranchAddress("metP4",&metP4);
    events->SetBranchAddress("metPtJESUP",&metPtJESUP);
    events->SetBranchAddress("metPtJESDOWN",&metPtJESDOWN);
    events->SetBranchAddress("metP4_GEN",&metP4_GEN);
    events->SetBranchAddress("metNoMu",&metNoMu);
    events->SetBranchAddress("pfMet_e3p0",&pfMet_e3p0);
    events->SetBranchAddress("trackMet",&trackMet);
    events->SetBranchAddress("caloMet_Pt",&caloMet_Pt);
    events->SetBranchAddress("caloMet_Phi",&caloMet_Phi);
    events->SetBranchAddress("caloMet_SumEt",&caloMet_SumEt);
    events->SetBranchAddress("genP4",&genP4);
    events->SetBranchAddress("genjetP4",&genjetP4);
    events->SetBranchAddress("genPdgId",&genPdgId);
    events->SetBranchAddress("puTrueInt",&puTrueInt);
    events->SetBranchAddress("mcWeight",&mcWeight);
    events->SetBranchAddress("pdfQscale",&pdfQscale);
    events->SetBranchAddress("pdfAlphaQED",&pdfAlphaQED);
    events->SetBranchAddress("pdfAlphaQCD",&pdfAlphaQCD);
    events->SetBranchAddress("pdfX1",&pdfX1);
    events->SetBranchAddress("pdfX2",&pdfX2);
    events->SetBranchAddress("pdfId1",&pdfId1);
    events->SetBranchAddress("pdfId2",&pdfId2);
    events->SetBranchAddress("pdfScalePdf",&pdfScalePdf);
    events->SetBranchAddress("photonP4",&photonP4);
    events->SetBranchAddress("photonIso",&photonIso);
    events->SetBranchAddress("photonSieie",&photonSieie);
    events->SetBranchAddress("photonTightId",&photonTightId);
    events->SetBranchAddress("photonMediumId",&photonMediumId);
    events->SetBranchAddress("photonLooseId",&photonLooseId);
    events->SetBranchAddress("photonChIso",&photonChIso);
    events->SetBranchAddress("photonChIsoRC",&photonChIsoRC);
    events->SetBranchAddress("photonNhIso",&photonNhIso);
    events->SetBranchAddress("photonNhIsoRC",&photonNhIsoRC);
    events->SetBranchAddress("photonPhoIso",&photonPhoIso);
    events->SetBranchAddress("photonPhoIsoRC",&photonPhoIsoRC);
    events->SetBranchAddress("photonPuIso",&photonPuIso);
    events->SetBranchAddress("photonPuIsoRC",&photonPuIsoRC);
    events->SetBranchAddress("tauP4",&tauP4);
    events->SetBranchAddress("tauId",&tauId);
    events->SetBranchAddress("tauQ",&tauQ);
    events->SetBranchAddress("tauM",&tauM);
    events->SetBranchAddress("tauIso",&tauIso);
    events->SetBranchAddress("triggerFired",&triggerFired);
    events->SetBranchAddress("triggerPrescale",&triggerPrescale);
    events->SetBranchAddress("triggerLeps",&triggerLeps);
    events->SetBranchAddress("triggerJets",&triggerJets);
    events->SetBranchAddress("triggerTaus",&triggerTaus);
    events->SetBranchAddress("triggerPhotons",&triggerPhotons);
    events->SetBranchAddress("npv",&npv);

    Long64_t nentries = events->GetEntries();
    unsigned int npairs=0;
    for (Long64_t i=0; i<nentries;i++) {
      events->GetEntry(i);
      unsigned int n_lep = lepPdgId->size();
      std::vector<TLorentzVector> p4_ele_tag_, p4_ele_passing_probe_, p4_ele_failing_probe_;
      std::vector<int> q_ele_tag_, q_ele_passing_probe_, q_ele_failing_probe_;
    
      // Object level loop
      if(n_lep != 0) { for(unsigned int il=0; il<n_lep; il++) {
        TLorentzVector *P4=(TLorentzVector*)lepP4->At(il);
        if( P4->Pt() >= 10.) { switch( abs( (*lepPdgId)[il] )) {
          case 11: // electron
            if(
              ((*lepSelBits)[il] & (0x1 << probe_id)) != 0 &&
              ((*lepSelBits)[il] & (0x1 << passing_probe_id)) == 0
            ) { // make probe selection but fail passing probe selection
              p4_ele_failing_probe_.push_back((*P4));
              if((*lepPdgId)[il] > 0) q_ele_failing_probe_.push_back(1);
              else q_ele_failing_probe_.push_back(-1);
            }
            if(((*lepSelBits)[il] & (0x1 << passing_probe_id)) != 0) { // pass tight ele ID i.e. passing probe or tag
              p4_ele_passing_probe_.push_back((*P4));
              if((*lepPdgId)[il] > 0) {
                q_ele_passing_probe_.push_back(1);
              } else {
                q_ele_passing_probe_.push_back(-1);
              }
            }
            if(
              ((*lepSelBits)[il] & (0x1 << tag_id)) != 0 
            ) { // pass tag ele ID, and trigger matching
              p4_ele_tag_.push_back((*P4));
              if((*lepPdgId)[il] > 0) {
                q_ele_tag_.push_back(1);
              } else {
                q_ele_tag_.push_back(-1);
              }
            }
            break;
          default:
            break;
        }}
      }}
  
      //demote some integers to unsigned because Kevin is a computer scientist
      out_runNum=runNum;
      out_lumiSec=lumiNum;
      out_evtNum=eventNum;
      out_npv=npv;
      npu=puTrueInt;
      
      scale1fb = h_npv_ratio->GetBinContent(h_npv_ratio->GetBin(npv));
      // associate pairs: electrons
      for(unsigned int iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
        // passing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_passing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            // if(verbose) printf("\t\tmade a PASSING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            npairs++;
            reweighted_tt_tree->Fill();
          }
        }
        // failing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_failing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            // if(verbose) printf("\t\tmade a FAILING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            npairs++;
            reweighted_tt_tree->Fill();
          }
        }

      }
    }
    if(verbose) printf("associated %d pairs\n",npairs);
    input_file->Close();
  }
  
  reweighted_tt_file->cd();
  reweighted_tt_tree->Write();
  reweighted_tt_file->Close();

  
  // make DY skim
  printf("making DY skim\n");
  while(getline(dy_ifs,input_file_name)) {
    if(verbose) printf("Opening file \"%s\"\n",input_file_name.c_str());
    TFile *input_file=TFile::Open(input_file_name.c_str(),"READ");
    TTree *events=(TTree*)input_file->Get("nero/events");
    // Set branch addresses.
    events->SetBranchAddress("isRealData",&isRealData);
    events->SetBranchAddress("runNum",&runNum);
    events->SetBranchAddress("lumiNum",&lumiNum);
    events->SetBranchAddress("eventNum",&eventNum);
    events->SetBranchAddress("rho",&rho);
    events->SetBranchAddress("jetP4",&jetP4);
    events->SetBranchAddress("jetRawPt",&jetRawPt);
    events->SetBranchAddress("jetBdiscr",&jetBdiscr);
    events->SetBranchAddress("jetBdiscrLegacy",&jetBdiscrLegacy);
    events->SetBranchAddress("jetPuId",&jetPuId);
    events->SetBranchAddress("jetUnc",&jetUnc);
    events->SetBranchAddress("jetQGL",&jetQGL);
    events->SetBranchAddress("jetFlavour",&jetFlavour);
    events->SetBranchAddress("jetMatchedPartonPdgId",&jetMatchedPartonPdgId);
    events->SetBranchAddress("jetMotherPdgId",&jetMotherPdgId);
    events->SetBranchAddress("jetGrMotherPdgId",&jetGrMotherPdgId);
    events->SetBranchAddress("jetMonojetId",&jetMonojetId);
    events->SetBranchAddress("jetMonojetIdLoose",&jetMonojetIdLoose);
    events->SetBranchAddress("jetQ",&jetQ);
    events->SetBranchAddress("jetQnoPU",&jetQnoPU);
    events->SetBranchAddress("lepP4",&lepP4);
    events->SetBranchAddress("lepPdgId",&lepPdgId);
    events->SetBranchAddress("lepIso",&lepIso);
    events->SetBranchAddress("lepSelBits",&lepSelBits);
    events->SetBranchAddress("lepPfPt",&lepPfPt);
    events->SetBranchAddress("lepChIso",&lepChIso);
    events->SetBranchAddress("lepNhIso",&lepNhIso);
    events->SetBranchAddress("lepPhoIso",&lepPhoIso);
    events->SetBranchAddress("lepPuIso",&lepPuIso);
    events->SetBranchAddress("metP4",&metP4);
    events->SetBranchAddress("metPtJESUP",&metPtJESUP);
    events->SetBranchAddress("metPtJESDOWN",&metPtJESDOWN);
    events->SetBranchAddress("metP4_GEN",&metP4_GEN);
    events->SetBranchAddress("metNoMu",&metNoMu);
    events->SetBranchAddress("pfMet_e3p0",&pfMet_e3p0);
    events->SetBranchAddress("trackMet",&trackMet);
    events->SetBranchAddress("caloMet_Pt",&caloMet_Pt);
    events->SetBranchAddress("caloMet_Phi",&caloMet_Phi);
    events->SetBranchAddress("caloMet_SumEt",&caloMet_SumEt);
    events->SetBranchAddress("genP4",&genP4);
    events->SetBranchAddress("genjetP4",&genjetP4);
    events->SetBranchAddress("genPdgId",&genPdgId);
    events->SetBranchAddress("puTrueInt",&puTrueInt);
    events->SetBranchAddress("mcWeight",&mcWeight);
    events->SetBranchAddress("pdfQscale",&pdfQscale);
    events->SetBranchAddress("pdfAlphaQED",&pdfAlphaQED);
    events->SetBranchAddress("pdfAlphaQCD",&pdfAlphaQCD);
    events->SetBranchAddress("pdfX1",&pdfX1);
    events->SetBranchAddress("pdfX2",&pdfX2);
    events->SetBranchAddress("pdfId1",&pdfId1);
    events->SetBranchAddress("pdfId2",&pdfId2);
    events->SetBranchAddress("pdfScalePdf",&pdfScalePdf);
    events->SetBranchAddress("photonP4",&photonP4);
    events->SetBranchAddress("photonIso",&photonIso);
    events->SetBranchAddress("photonSieie",&photonSieie);
    events->SetBranchAddress("photonTightId",&photonTightId);
    events->SetBranchAddress("photonMediumId",&photonMediumId);
    events->SetBranchAddress("photonLooseId",&photonLooseId);
    events->SetBranchAddress("photonChIso",&photonChIso);
    events->SetBranchAddress("photonChIsoRC",&photonChIsoRC);
    events->SetBranchAddress("photonNhIso",&photonNhIso);
    events->SetBranchAddress("photonNhIsoRC",&photonNhIsoRC);
    events->SetBranchAddress("photonPhoIso",&photonPhoIso);
    events->SetBranchAddress("photonPhoIsoRC",&photonPhoIsoRC);
    events->SetBranchAddress("photonPuIso",&photonPuIso);
    events->SetBranchAddress("photonPuIsoRC",&photonPuIsoRC);
    events->SetBranchAddress("tauP4",&tauP4);
    events->SetBranchAddress("tauId",&tauId);
    events->SetBranchAddress("tauQ",&tauQ);
    events->SetBranchAddress("tauM",&tauM);
    events->SetBranchAddress("tauIso",&tauIso);
    events->SetBranchAddress("triggerFired",&triggerFired);
    events->SetBranchAddress("triggerPrescale",&triggerPrescale);
    events->SetBranchAddress("triggerLeps",&triggerLeps);
    events->SetBranchAddress("triggerJets",&triggerJets);
    events->SetBranchAddress("triggerTaus",&triggerTaus);
    events->SetBranchAddress("triggerPhotons",&triggerPhotons);
    events->SetBranchAddress("npv",&npv);

    Long64_t nentries = events->GetEntries();
    unsigned int npairs=0;
    for (Long64_t i=0; i<nentries;i++) {
      events->GetEntry(i);
      unsigned int n_lep = lepPdgId->size();
      std::vector<TLorentzVector> p4_ele_tag_, p4_ele_passing_probe_, p4_ele_failing_probe_;
      std::vector<int> q_ele_tag_, q_ele_passing_probe_, q_ele_failing_probe_;
    
      // Object level loop
      if(n_lep != 0) { for(unsigned int il=0; il<n_lep; il++) {
        TLorentzVector *P4=(TLorentzVector*)lepP4->At(il);
        if( P4->Pt() >= 10.) { switch( abs( (*lepPdgId)[il] )) {
          case 11: // electron
            if(
              ((*lepSelBits)[il] & (0x1 << probe_id)) != 0 &&
              ((*lepSelBits)[il] & (0x1 << passing_probe_id)) == 0
            ) { // make probe selection but fail passing probe selection
              p4_ele_failing_probe_.push_back((*P4));
              if((*lepPdgId)[il] > 0) q_ele_failing_probe_.push_back(1);
              else q_ele_failing_probe_.push_back(-1);
            }
            if(((*lepSelBits)[il] & (0x1 << passing_probe_id)) != 0) { // pass tight ele ID i.e. passing probe or tag
              p4_ele_passing_probe_.push_back((*P4));
              if((*lepPdgId)[il] > 0) {
                q_ele_passing_probe_.push_back(1);
              } else {
                q_ele_passing_probe_.push_back(-1);
              }
            }
            if(
              ((*lepSelBits)[il] & (0x1 << tag_id)) != 0
            ) { // pass tag ele ID, and trigger matching
              p4_ele_tag_.push_back((*P4));
              if((*lepPdgId)[il] > 0) {
                q_ele_tag_.push_back(1);
              } else {
                q_ele_tag_.push_back(-1);
              }
            }
            break;
          default:
            break;
        }}
      }}
  
      //demote some integers to unsigned because Kevin is a computer scientist
      out_runNum=runNum;
      out_lumiSec=lumiNum;
      out_evtNum=eventNum;
      out_npv=npv;
      npu=puTrueInt;
      
      scale1fb = h_npv_ratio->GetBinContent(h_npv_ratio->GetBin(npv));
      // associate pairs: electrons
      for(unsigned int iTag=0; iTag < p4_ele_tag_.size(); iTag++) {
        // passing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_passing_probe_.size(); iProbe++) {
          pass = 1;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_passing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_passing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a PASSING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            reweighted_dy_tree->Fill();
            npairs++;
          }
        }
        // failing probes for electrons
        for(unsigned int iProbe=0; iProbe < p4_ele_failing_probe_.size(); iProbe++) {
          pass = 0;
          *p4_tag     = p4_ele_tag_[iTag];
          *p4_probe   = p4_ele_failing_probe_[iProbe];
          qtag    = q_ele_tag_[iTag];
          qprobe  = q_ele_failing_probe_[iProbe];
          // make sure they aren't the same particle with tiny delta R veto
          if(!( pow((p4_tag->Pt() - p4_probe->Pt()), 2) + pow((p4_tag->Eta() - p4_probe->Eta()), 2) < .000001) ) {
            TLorentzVector systemP4 = (*p4_tag) + (*p4_probe);
            mass = systemP4.M();
            if(verbose) printf("\t\tmade a FAILING electron pair! pTs %f, %f; system mass %f, total charge %d e\n", p4_tag->Pt(), p4_probe->Pt(), mass, qtag+qprobe);
            reweighted_dy_tree->Fill();
            npairs++;
          }
        }

      }
    }
    if(verbose) printf("associated %d pairs\n",npairs);
    input_file->Close();
  }
  
  reweighted_dy_file->cd();
  reweighted_dy_tree->Write();
  reweighted_dy_file->Close();

}

