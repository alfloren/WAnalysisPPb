#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include <TLorentzVector.h>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include <vector>
#include <string>

const Int_t MAXPARTICLE = 100000;
const Int_t MAXJET = 100;
//
// DiJet ana Event Data Tree definition
//
class TreePFCandEventData
{
 public:
  // ===== Class Methods =====
  void SetDefaults();
  TreePFCandEventData();
  void SetTree(TTree * t) { tree_=t; }
  void SetBranches();
  void Clear();
  bool doJets;
  bool doMC;
   Int_t centBin, hiNtracks; 
  Double_t c;
  Int_t nEvent, nRun, nLumi;
  Float_t recoPfMET_, recoPfMETPhi_, recoPfMETsumEt_, recoPfMETmEtSig_, recoPfMETSig_;
  Float_t VertexZ;

 Float_t hiHF, hiHFplus, hiHFminus, hiHFplusEta4, hiHFminusEta4;

  Float_t                 jdphi_;
  // -- particle info --
  bool                    isTrackerMuon_[MAXPARTICLE], isGlobalMuon_[MAXPARTICLE], pass12Filter_[MAXPARTICLE], pass7Filter_[MAXPARTICLE], pass3Filter_[MAXPARTICLE], pass12Path_[MAXPARTICLE], pass7Path_[MAXPARTICLE], pass3Path_[MAXPARTICLE], ;
  Int_t                   nPFpart_, nGENpart_, njets_, nGENNeu_, nGENW_, nTRACKpart_,nMuonpart_, nPFmuon_, nGENMu_, nGMpart_, nGENWout_;
  Int_t                   pfId_[MAXPARTICLE], genPDGId_[MAXPARTICLE], traQual_[MAXPARTICLE], genMotherId_[20], genTauCharge_; 
  Int_t                   genMotherStatus_[MAXPARTICLE], traAlgo_[MAXPARTICLE], pfTrackhits_[MAXPARTICLE], trahits_[MAXPARTICLE], genWstatus_, genWcharge_, genNbMother_[MAXPARTICLE], genStatus_[MAXPARTICLE];
  Int_t                   pfCharge_[MAXPARTICLE], genCharge_[MAXPARTICLE], muCharge_[MAXPARTICLE],genMuCharge_[2], pfTrackerMuon_[MAXPARTICLE], pfdxy_[MAXPARTICLE],  pfchi2_[MAXPARTICLE], muTrackHits_[MAXPARTICLE], muPixHits_[MAXPARTICLE], muchi2_[MAXPARTICLE], muSeg_[MAXPARTICLE], nMuValHits_[MAXPARTICLE], chi2Track_[MAXPARTICLE], StandAloneMuon_[MAXPARTICLE];
  Float_t                 pfPt_[MAXPARTICLE], genPt_[MAXPARTICLE],  jetPt_[MAXPARTICLE], pfGloMuon_[MAXPARTICLE], GMpatPt_[MAXPARTICLE], genTauPt_, genTauEta_, genTauPhi_, genTaupx_, genTaupy_;
  Float_t                 pfEta_[MAXPARTICLE], genEta_[MAXPARTICLE],  jetEta_[MAXPARTICLE], genWmass_;
  Float_t                 jetMass_[MAXJET], jetY_[MAXJET], jetPU_[MAXJET], rawPt_[MAXJET], genMupx_[MAXPARTICLE], genMupy_[MAXPARTICLE];
   Float_t                 pfPhi_[MAXPARTICLE], genPhi_[MAXPARTICLE],  jetPhi_[MAXPARTICLE], genWPt_, genWEta_, genWPhi_, genWrapidity_;
  Float_t                 pfTheta_[MAXPARTICLE], genTheta_[MAXPARTICLE],  jetTheta_[MAXPARTICLE]; 
  Float_t                 traPhi_[MAXPARTICLE], traEta_[MAXPARTICLE],traPt_[MAXPARTICLE], muPt_[MAXPARTICLE], muPhi_[MAXPARTICLE], muEta_[MAXPARTICLE], mutrackPt_[MAXPARTICLE], muPx_[MAXPARTICLE], muPy_[MAXPARTICLE], muPz_[MAXPARTICLE], muTrackIso_[MAXPARTICLE], muCaloIso_[MAXPARTICLE], muEcalIso_[MAXPARTICLE], muHcalIso_[MAXPARTICLE];
  Float_t                 genNeuEt_, genNeuEta_, genNeuPhi_, genNeuPt_, genMuPt_[2], genMuEta_[2], genMuPhi_[2], pfMuonpx_[MAXPARTICLE],pfMuonpy_[MAXPARTICLE],pfMuonpz_[MAXPARTICLE],pfMuone_[MAXPARTICLE], pfChargepx_[MAXPARTICLE],pfChargepy_[MAXPARTICLE],pfChargepz_[MAXPARTICLE] ; 
  Float_t                 pfEt_[MAXPARTICLE], genEt_[MAXPARTICLE],jetEt_[MAXPARTICLE], pfTrack_[MAXPARTICLE], pfTrackdijet_[MAXPARTICLE] ;
  Int_t                   jetMatchIndex_[MAXPARTICLE];
Float_t  mudxy_[MAXPARTICLE];
Float_t stdTrack_[MAXPARTICLE];  

  
 private:
  TTree*                 tree_;

};

class PFCandAnalyzer : public edm::EDAnalyzer {
  public:
    explicit PFCandAnalyzer(const edm::ParameterSet&);
    ~PFCandAnalyzer();

    // class methods


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

   


     bool findW(const reco::Candidate *gen);
     bool findWbis(const reco::Candidate *gen);
      bool findWTris(const reco::Candidate *gen);
    // ----------member data ---------------------------
    edm::Service<TFileService> fs;

    // === Ana setup ===

    // Event Info
    edm::InputTag pfCandidateLabel_;
    edm::InputTag genLabel_;
    edm::InputTag jetLabel_;
    edm::InputTag trackLabel_;
    edm::InputTag muonLabel_;
    edm::InputTag recoPfMETlabel_;
    edm::InputTag _thePVs;
    
    TFile *fhisto;
    TH1F *beamspot; 
    TTree	  *pfTree_;
    TreePFCandEventData pfEvt_;

    // cuts
    Double_t        pfPtMin_;
    Double_t        jetPtMin_;
    Double_t        genPtMin_;
    Double_t        dxyCut_;
   Double_t        normalizedChi2Cut_;
   
   Double_t        ptThrForZ1_;
   Double_t         ptThrForZ2_;
   bool            isRelativeIso_;
   bool            isCombinedIso_;
   double           isoCut03_;
   
   Int_t trackerHitsCut_;
  Int_t pixelHitsCut_;
  bool isAlsoTrackerMuon_;
   // debug
    Int_t	  verbosity_;
    bool   useQuality_;
    bool   doJets_;
    bool   doMC_;
    bool   skipCharged_;
    
    std::string qualityString_;
    //centrality
 CentralityProvider * centrality_;
};

