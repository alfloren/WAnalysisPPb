// -*- C++ -*-
//
// Package:    HiPFCandAnalyzer
// Class:      HiPFCandAnalyzer
// 
/**\class HiPFCandAnalyzer HiPFCandAnalyzer.cc ana/HiPFCandAnalyzer/src/HiPFCandAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Matt, Nguyen 
//         Created:   Oct  10 2010
// 
//
//


// system include files
#include <memory>

// stl
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ana
#include "alfloren/pfMetAnalysis/interface/PFCandAnalyzer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include  "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

//
// constructors and destructor
//
PFCandAnalyzer::PFCandAnalyzer(const edm::ParameterSet& iConfig)
{
  // Event source
  // Event Info
  pfCandidateLabel_ = iConfig.getParameter<edm::InputTag>("pfCandidateLabel");
  genLabel_ = iConfig.getParameter<edm::InputTag>("genLabel");
  jetLabel_ = iConfig.getParameter<edm::InputTag>("jetLabel");
  trackLabel_ = iConfig.getUntrackedParameter<edm::InputTag> ("TrackLabel", edm::InputTag("generalTracks"));
  muonLabel_ = iConfig.getUntrackedParameter<edm::InputTag> ("MuonLabel",edm::InputTag("goodPatMuons"));
  // muonLabel_ = iConfig.getUntrackedParameter<edm::InputTag> ("MuonLabel",edm::InputTag("remuons"));
  recoPfMETlabel_  = iConfig.getParameter<edm::InputTag>("recoPfMETLabel");
  pfPtMin_ = iConfig.getParameter<double>("pfPtMin");
  genPtMin_ = iConfig.getParameter<double>("genPtMin");
  jetPtMin_ = iConfig.getParameter<double>("jetPtMin");
  _thePVs = iConfig.getParameter<edm::InputTag>("VtxLabel");

  // debug
  verbosity_ = iConfig.getUntrackedParameter<int>("verbosity", 0);
  useQuality_ = iConfig.getUntrackedParameter<bool>("useQuality",false);
  doJets_ = iConfig.getUntrackedParameter<bool>("doJets",0);
  doMC_ = iConfig.getUntrackedParameter<bool>("doMC",0);
  skipCharged_ = iConfig.getUntrackedParameter<bool>("skipCharged",0);
  qualityString_ = iConfig.getUntrackedParameter<std::string>("qualityString","highPurity");

  //muon cuts
  dxyCut_ = (iConfig.getUntrackedParameter<double>("DxyCut", 0.02)); // rejeter les muons cosmiques
  
  normalizedChi2Cut_ = (iConfig.getUntrackedParameter<double>("NormalizedChi2Cut", 10.)); // rejeter les muons venant du decay des hadrons
  
  isRelativeIso_ = iConfig.getUntrackedParameter<bool>("IsRelativeIso", true);
  isCombinedIso_ = iConfig.getUntrackedParameter<bool>("IsCombinedIso", false);
  isoCut03_ = iConfig.getUntrackedParameter<double>("IsoCut03", 0.15);
  trackerHitsCut_ = iConfig.getUntrackedParameter<int>("TrackerHitsCut", 11);
  pixelHitsCut_ = iConfig.getUntrackedParameter<int>("PixelHitsCut", 1);
  isAlsoTrackerMuon_ = iConfig.getUntrackedParameter<bool>("IsAlsoTrackerMuon", true);
  // Z rejection
  ptThrForZ1_ = (iConfig.getUntrackedParameter<double>("PtThrForZ1", 20.));
  ptThrForZ2_ = (iConfig.getUntrackedParameter<double>("PtThrForZ2", 10.));


  centrality_= 0;
 
  fhisto = new TFile("myhisto.root","RECREATE");
  fhisto->cd();
  beamspot = new TH1F("beamspot","",20,0,20);
  beamspot->SetDirectory(fhisto);

}


PFCandAnalyzer::~PFCandAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------


 
 void
PFCandAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

 
//centrality

   
 pfEvt_.Clear();

 

if(!centrality_) centrality_ = new CentralityProvider(iSetup);
   centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
   //double c = centrality_->centralityValue();
  const reco::Centrality* centrality =  centrality_->raw();
  
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  // Fill PF info
  // cout<<"beamspot ok"<<endl;
  
 if (!iEvent.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle)) {
   cout<<" No beam spot found !!!"<<endl;
  }
 

   edm::Handle<reco::VertexCollection> privtxs;
   iEvent.getByLabel(_thePVs, privtxs);

   reco::VertexCollection::const_iterator privtx;
   math::XYZPoint RefVtx;
   float zVtx;
   float nPV;

     if ( privtxs->begin() != privtxs->end() ) {
    privtx=privtxs->begin();
    RefVtx = privtx->position();
  } else {
    RefVtx.SetXYZ(0.,0.,0.);
  }

    zVtx = RefVtx.Z();

  pfEvt_.hiHF = centrality->EtHFtowerSum();
  pfEvt_.hiHFplus = centrality->EtHFtowerSumPlus();
   pfEvt_.hiHFminus = centrality->EtHFtowerSumMinus();
  pfEvt_.hiHFplusEta4 = centrality->EtHFtruncatedPlus();
  pfEvt_.hiHFminusEta4 = centrality->EtHFtruncatedMinus();
  pfEvt_.hiNtracks = centrality->Ntracks();
 pfEvt_.centBin = centrality_->getBin();
 pfEvt_.nEvent = iEvent.id().event();
 pfEvt_.nRun = iEvent.id().run();
 pfEvt_.nLumi = iEvent.luminosityBlock();
 pfEvt_.VertexZ = zVtx;


  // Fill Jet info
 /* edm::Handle<pat::JetCollection> jets;
 iEvent.getByLabel(jetLabel_,jets);  
 const pat::JetCollection *jetColl = &(*jets);
 */
 /* if(doJets_){
     
     for(unsigned ijet=0;ijet<jetColl->size(); ijet++) {
       const pat::Jet jet = jetColl->at(ijet);
       
       double pt =  jet.pt();
       
       if(pt>jetPtMin_){
	 pfEvt_.rawPt_[pfEvt_.njets_]=jet.correctedJet("Uncorrected").pt();    
	 pfEvt_.jetPt_[pfEvt_.njets_] = pt;      
	 pfEvt_.jetEta_[pfEvt_.njets_] = jet.eta();      
	 pfEvt_.jetPhi_[pfEvt_.njets_] = jet.phi();      
	 pfEvt_.jetMass_[pfEvt_.njets_] = jet.mass();   
	 pfEvt_.jetY_[pfEvt_.njets_] = jet.eta();
	 pfEvt_.jetPU_[pfEvt_.njets_] = jet.pileup();   

	 pfEvt_.njets_++;
       }
     }
     }	
 */

 //Fill PF candidates
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByLabel(pfCandidateLabel_,pfCandidates);  
  const reco::PFCandidateCollection *pfCandidateColl = &(*pfCandidates);

  
  for(unsigned icand=0;icand<pfCandidateColl->size(); icand++) {
    const reco::PFCandidate pfCandidate = pfCandidateColl->at(icand);      





      double pt =  pfCandidate.pt();
      if(pt<pfPtMin_) continue;

      int id = pfCandidate.particleId();
      if(skipCharged_ && (abs(id) == 1 || abs(id) == 3)) continue;

      bool matched2Jet = false;

      /*  if(doJets_){
	
	pfEvt_.jetMatchIndex_[pfEvt_.nPFpart_] = -1;

	for(unsigned ijet=0;ijet<jetColl->size(); ijet++) {
	  const pat::Jet jet = jetColl->at(ijet);
	  	  
	  if(jet.pt()>jetPtMin_){

	    std::vector<reco::PFCandidatePtr> pfConstituents = jet.getPFConstituents();
	    
	    for( std::vector<reco::PFCandidatePtr>::const_iterator ibegin=pfConstituents.begin(), iend=pfConstituents.end(), iconstituent=ibegin; iconstituent!=iend; ++iconstituent){
	      
	      reco::PFCandidatePtr candptr(pfCandidates, icand);
		edm::Ptr<reco::PFCandidate> pfBackRef ( *iconstituent );
		
		// couldn't figure out the matching by ref, so just do it like this:
		if(pfBackRef->pt() == pfCandidate.pt() && pfBackRef->eta()== pfCandidate.eta() && pfBackRef->particleId()== pfCandidate.particleId() ){
		  
		  pfEvt_.jetMatchIndex_[pfEvt_.nPFpart_] = ijet;
		  matched2Jet =true;
		  break;
		}
	    }

	  }
	  if(matched2Jet== true) break;
	}
	}*/
      

      pfEvt_.pfId_[pfEvt_.nPFpart_] = id;      
      pfEvt_.pfPt_[pfEvt_.nPFpart_] = pt;      
      pfEvt_.pfEta_[pfEvt_.nPFpart_] = pfCandidate.eta();      
      pfEvt_.pfPhi_[pfEvt_.nPFpart_] = pfCandidate.phi(); 
      pfEvt_.pfTheta_[pfEvt_.nPFpart_] = pfCandidate.theta(); 
      pfEvt_.pfEt_[pfEvt_.nPFpart_] = pfCandidate.et();  
      pfEvt_.pfCharge_[pfEvt_.nPFpart_] = pfCandidate.charge(); 
      pfEvt_.pfTrack_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfTrackerMuon_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfdxy_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfchi2_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfMuonpx_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfMuonpy_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfMuonpz_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfChargepy_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfChargepx_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfChargepz_[pfEvt_.nPFpart_] = 0;
      pfEvt_.pfTrackdijet_[pfEvt_.nPFpart_] = 0;
      if(abs(id) == 3){
      
	pfEvt_.pfMuonpx_[pfEvt_.nPFpart_] = pfCandidate.px();
        pfEvt_.pfMuonpy_[pfEvt_.nPFpart_] = pfCandidate.py();
        pfEvt_.pfMuonpz_[pfEvt_.nPFpart_] = pfCandidate.pz();
	//pfEvt_.pfMuone_[pfEvt_.nPFmuon_] = pfCandidate.e();
	//cout<<pfCandidate.px()<<endl;
      const reco::MuonRef muonRef = pfCandidate.muonRef();  
      //const reco::Muon& muon = *muonRef;
      if ( muonRef->isTrackerMuon() ){
      reco::TrackRef trackRef = muonRef->track();
      
      double dxy = trackRef->dxy(beamSpotHandle->position());
      
      double normalizedChi2 = trackRef->normalizedChi2();
       
      // const reco::Track& track = *trackRef;
      pfEvt_.pfTrack_[pfEvt_.nPFpart_] = trackRef->pt();
      pfEvt_.pfTrackhits_[pfEvt_.nPFpart_] = trackRef->numberOfValidHits();
      //if(trackRef->numberOfValidHits()<=12) cout<<"mauvais hits"<<trackRef->numberOfValidHits()<<endl;
      pfEvt_.pfTrackerMuon_[pfEvt_.nPFpart_] = 1;
      if (fabs(dxy)<dxyCut_) pfEvt_.pfdxy_[pfEvt_.nPFpart_] = 1; 
      if(fabs(normalizedChi2)<normalizedChi2Cut_) pfEvt_.pfchi2_[pfEvt_.nPFpart_] = 1; 

     if ( muonRef->isGlobalMuon() ){
        reco::TrackRef globalmuon = muonRef->globalTrack();
        pfEvt_.pfGloMuon_[pfEvt_.nGMpart_] = globalmuon->pt();
       pfEvt_.nGMpart_++;
      
   }   
 }


  
      pfEvt_.nPFmuon_++;
      }
      

      if(abs(id) == 1){
      pfEvt_.pfChargepx_[pfEvt_.nPFpart_] = pfCandidate.px();
      pfEvt_.pfChargepy_[pfEvt_.nPFpart_] = pfCandidate.py();
      pfEvt_.pfChargepz_[pfEvt_.nPFpart_] = pfCandidate.pz();
      const reco::TrackRef trackRef = pfCandidate.trackRef(); 
      pfEvt_.pfTrackdijet_[pfEvt_.nPFpart_] = trackRef->pt();
      
      
       }
      pfEvt_.nPFpart_++;

      
      
  }
	
  

  // Fill GEN info
    if(doMC_){
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genLabel_,genParticles);     
  const reco::GenParticleCollection* genColl= &(*genParticles);
 
  for(unsigned igen=0;igen<genColl->size(); igen++) {    

    const reco::GenParticle gen = genColl->at(igen);    
   
   

    
 
    if(abs(gen.pdgId())==16 && gen.status()==1){

 
	//cout<<" I have a neutrino "<<endl;
	const Candidate *neutrino = &gen;
	
		bool isFromW = findW(neutrino);

        //cout<<"++++++++++"<<isFromW<<endl;
	if (isFromW){
	  //std::cout << " Et neutrino= " <<gen.et()<<std::endl;
	  pfEvt_.genNeuEt_ = gen.et();
          pfEvt_.genNeuPt_ = gen.pt();   
          pfEvt_.genNeuEta_ =gen.eta();
	  pfEvt_.genNeuPhi_ =gen.phi();
          	pfEvt_.nGENNeu_ ++ ;
	  }	
      }

    if(abs(gen.pdgId()) == 13 && gen.status()==1){

	
	const Candidate *muon = &gen;

        bool isFromWbis = findWbis(muon);

        if (isFromWbis){
	 
	  pfEvt_.genMuPt_[pfEvt_.nGENMu_] = gen.pt();  
          pfEvt_.genMuEta_[pfEvt_.nGENMu_] =gen.eta();
	  pfEvt_.genMuPhi_[pfEvt_.nGENMu_] =gen.phi();
          pfEvt_.genMuCharge_[pfEvt_.nGENMu_] =gen.charge();
          pfEvt_.genMupx_[pfEvt_.nGENMu_] =gen.px();
          pfEvt_.genMupy_[pfEvt_.nGENMu_] =gen.py();
      pfEvt_.nGENMu_++;
	  }	

    }
   
    if(abs(gen.pdgId()) == 15 && gen.status()==2){

	
	const Candidate *tau = &gen;

        bool isFromWTris = findWTris(tau);

        if (isFromWTris){
	 
	  pfEvt_.genTauPt_ = gen.pt();  
          pfEvt_.genTauEta_ =gen.eta();
	  pfEvt_.genTauPhi_ =gen.phi();
          pfEvt_.genTauCharge_ =gen.charge();
          pfEvt_.genTaupx_ =gen.px();
          pfEvt_.genTaupy_ =gen.py();
    
	  }	

	  }

    if( gen.pt()> genPtMin_){      
    
      pfEvt_.genPDGId_[pfEvt_.nGENpart_] = gen.pdgId();
      pfEvt_.genStatus_[pfEvt_.nGENpart_] = gen.status();        
      pfEvt_.genPt_[pfEvt_.nGENpart_] = gen.pt();      
      pfEvt_.genEta_[pfEvt_.nGENpart_] = gen.eta();      
      pfEvt_.genPhi_[pfEvt_.nGENpart_] = gen.phi();
      pfEvt_.genEt_[pfEvt_.nGENpart_] = gen.et();  
      pfEvt_.genCharge_[pfEvt_.nGENpart_] = gen.charge();  
      pfEvt_.genNbMother_[pfEvt_.nGENpart_] = gen.numberOfMothers(); 

     
      Int_t nMother = gen.numberOfMothers();
      for(Int_t y=0; y<nMother; y++){    
       
      pfEvt_.genMotherId_[y] = gen.mother()->pdgId();
      
      }
       pfEvt_.nGENpart_++;
    }
    // pdgid= 23 is Z, pdgid =24 is w
    if(abs(gen.pdgId())==23){
          pfEvt_.genWPt_ = gen.pt();   
          pfEvt_.genWEta_ =gen.eta();
	  pfEvt_.genWPhi_ =gen.phi();    
          pfEvt_.genWstatus_ =gen.status();
          pfEvt_.genWcharge_ =gen.charge();
          pfEvt_.genWrapidity_ =gen.rapidity();
          pfEvt_.genWmass_ = gen.mass();
         
       }    
    
  }
  }
  
  //std::cout << " nbre de gen = " <<pfEvt_.nGENpart_<<std::endl;


  // Fill tracks

  edm::Handle<reco::TrackCollection> trackCollection;
 iEvent.getByLabel(trackLabel_, trackCollection);
 //cout << " TRACKS Size " << trackCollection->size() << endl;

const reco::TrackCollection* trackColl= &(*trackCollection);
 
  for(unsigned itrack=0;itrack<trackColl->size(); itrack++) { 

      const reco::Track tra = trackColl->at(itrack); 
       
      pfEvt_.traQual_[pfEvt_.nTRACKpart_] = 0;
      if(tra.quality(reco::TrackBase::qualityByName(qualityString_))) pfEvt_.traQual_[pfEvt_.nTRACKpart_]=1;
      
      //if(useQuality_ && pfEvt_.traQual_[pfEvt_.nTRACKpart_] != 1) continue;  
  
      
      pfEvt_.traPt_[pfEvt_.nTRACKpart_] = tra.pt();      
      pfEvt_.traEta_[pfEvt_.nTRACKpart_] = tra.eta();      
      pfEvt_.traPhi_[pfEvt_.nTRACKpart_] = tra.phi();
      pfEvt_.traAlgo_[pfEvt_.nTRACKpart_] = tra.algo();
      pfEvt_.trahits_[pfEvt_.nTRACKpart_] = tra.numberOfValidHits();
      // if(pfEvt_.traQual_[pfEvt_.nTRACKpart_]==1 && tra.numberOfValidHits()<=12) cout<<"number of hits pas bon ="<<tra.numberOfValidHits()<<endl;
      pfEvt_.nTRACKpart_++;

   }


  //MET Object

  Handle<reco::PFMETCollection>   recoPfMETHandle;
    iEvent.getByLabel(recoPfMETlabel_, recoPfMETHandle);

    const reco::PFMET& pfmet = recoPfMETHandle->at(0);

    if (recoPfMETHandle.isValid()) {
      pfEvt_.recoPfMET_ = pfmet.et();
      pfEvt_.recoPfMETPhi_ = pfmet.phi();
      pfEvt_.recoPfMETsumEt_  = pfmet.sumEt();
      pfEvt_.recoPfMETmEtSig_ = pfmet.mEtSig();
      pfEvt_.recoPfMETSig_    = pfmet.significance();
      } 
   
   
  //fill muon info
  //edm::Handle<reco::MuonCollection> muonCollection;//commenter pour le pat
  Handle< View<pat::Muon> > muonCollection; 
   iEvent.getByLabel(muonLabel_, muonCollection);//commenter pour le pat
 // cout << " muon Size " << muonCollection->size() << endl;

 //const reco::MuonCollection* muonColl= &(*muonCollection);//commenter pour le pat
 
  
   for(unsigned imuon=0;imuon<muonCollection->size(); imuon++) { 
     // const reco::Muon mu = muonColl->at(imuon);//commenter pour le pat
      
     
     const pat::Muon& mu = muonCollection->at(imuon);
     
    

     pfEvt_.muPt_[pfEvt_.nMuonpart_] = mu.pt(); 
      pfEvt_.muPx_[pfEvt_.nMuonpart_] = mu.px();
      pfEvt_.muPy_[pfEvt_.nMuonpart_] = mu.py();
      pfEvt_.muPz_[pfEvt_.nMuonpart_] = mu.pz();     
      pfEvt_.muEta_[pfEvt_.nMuonpart_] = mu.eta();      
      pfEvt_.muPhi_[pfEvt_.nMuonpart_] = mu.phi();
      pfEvt_.muCharge_[pfEvt_.nMuonpart_] = mu.charge();
      pfEvt_.muTrackIso_[pfEvt_.nMuonpart_] = mu.trackIso();
      pfEvt_.muCaloIso_[pfEvt_.nMuonpart_] = mu.caloIso();
      pfEvt_.muEcalIso_[pfEvt_.nMuonpart_] = mu.ecalIso();
      pfEvt_.muHcalIso_[pfEvt_.nMuonpart_] = mu.hcalIso();

       

      const pat::TriggerObjectStandAloneCollection muHLTMatchesFilter12 = mu.triggerObjectMatchesByFilter("hltL3fL2sMu12L3Filtered12");
      const pat::TriggerObjectStandAloneCollection muHLTMatchesPath12 = mu.triggerObjectMatchesByPath("HLT_PAMu12_v1",true,false); 
      const pat::TriggerObjectStandAloneCollection muHLTMatchesFilter7 = mu.triggerObjectMatchesByFilter("hltL3fL2sMu7L3Filtered7");
      const pat::TriggerObjectStandAloneCollection muHLTMatchesPath7 = mu.triggerObjectMatchesByPath("HLT_PAMu7_v1",true,false);
       const pat::TriggerObjectStandAloneCollection muHLTMatchesFilter3 = mu.triggerObjectMatchesByFilter("hltL3fL2sMu3L3Filtered3");
       const pat::TriggerObjectStandAloneCollection muHLTMatchesPath3 = mu.triggerObjectMatchesByPath("HLT_PAMu3_v1",true,false);


        pfEvt_.pass12Filter_[pfEvt_.nMuonpart_] = false;
        pfEvt_.pass7Filter_[pfEvt_.nMuonpart_] = false;
        pfEvt_.pass3Filter_[pfEvt_.nMuonpart_] = false;
        pfEvt_.pass12Path_[pfEvt_.nMuonpart_] = false;
        pfEvt_.pass7Path_[pfEvt_.nMuonpart_] = false;
        pfEvt_.pass3Path_[pfEvt_.nMuonpart_] = false;



	pfEvt_.pass12Filter_[pfEvt_.nMuonpart_] = (muHLTMatchesFilter12.size() >0);
        pfEvt_.pass7Filter_[pfEvt_.nMuonpart_] = (muHLTMatchesFilter7.size() >0);
        pfEvt_.pass3Filter_[pfEvt_.nMuonpart_] = (muHLTMatchesFilter3.size() >0);
        pfEvt_.pass12Path_[pfEvt_.nMuonpart_] = (muHLTMatchesPath12.size() >0);
        pfEvt_.pass7Path_[pfEvt_.nMuonpart_] = (muHLTMatchesPath7.size() >0);
        pfEvt_.pass3Path_[pfEvt_.nMuonpart_] = (muHLTMatchesPath3.size() >0);
     
       


      // cout<<"ok"<<endl;

     pfEvt_.muTrackHits_[pfEvt_.nMuonpart_]=0;
     pfEvt_.muPixHits_[pfEvt_.nMuonpart_]=0;
     pfEvt_.mudxy_[pfEvt_.nMuonpart_]=10.;
     pfEvt_.chi2Track_[pfEvt_.nMuonpart_] = 20;
     pfEvt_.muchi2_[pfEvt_.nMuonpart_] =20;
     pfEvt_.muSeg_[pfEvt_.nMuonpart_] = 0;
     pfEvt_.nMuValHits_[pfEvt_.nMuonpart_] =0;

     const reco::TrackRef mutrack = mu.innerTrack();
     
     const reco::TrackRef gm = mu.globalTrack();

     const reco::TrackRef standAlone = mu.standAloneMuon();

    const reco::TrackRef standAlonetrack = mu.outerTrack();

     Int_t trackerHits=0;
     Int_t pixelHits =0;
   
     
  

    
     
     pfEvt_.StandAloneMuon_[pfEvt_.nMuonpart_] = 0; 
     if(standAlone.isNonnull()){
      pfEvt_.StandAloneMuon_[pfEvt_.nMuonpart_] = 1;
      pfEvt_.stdTrack_[pfEvt_.nMuonpart_] = standAlonetrack->pt();
      }
         pfEvt_.isGlobalMuon_[pfEvt_.nMuonpart_] = mu.isGlobalMuon();
         pfEvt_.isTrackerMuon_[pfEvt_.nMuonpart_] = mu.isTrackerMuon();

     if (gm.isNonnull()){
      pfEvt_.muchi2_[pfEvt_.nMuonpart_] = gm->normalizedChi2();  
      pfEvt_.muSeg_[pfEvt_.nMuonpart_] = mu.numberOfMatches();
       pfEvt_.GMpatPt_[pfEvt_.nMuonpart_] =gm->pt();   
       pfEvt_.nMuValHits_[pfEvt_.nMuonpart_] = gm->hitPattern().numberOfValidMuonHits();
      }


    
    if (mutrack.isNonnull()){
     

      // if (fabs(mutrack->dxy(beamSpotHandle->position()))<dxyCut_) pfEvt_.mudxy_[pfEvt_.nMuonpart_] = 1; 
    
     trackerHits = mutrack->hitPattern().numberOfValidTrackerHits();
     pixelHits = mutrack->hitPattern().numberOfValidPixelHits(); 
    
    pfEvt_.muTrackHits_[pfEvt_.nMuonpart_]=trackerHits;
    pfEvt_.muPixHits_[pfEvt_.nMuonpart_]=pixelHits;
    pfEvt_.mudxy_[pfEvt_.nMuonpart_]=mutrack->dxy(beamSpotHandle->position());
    pfEvt_.mutrackPt_[pfEvt_.nMuonpart_] = mutrack->pt();
    pfEvt_.chi2Track_[pfEvt_.nMuonpart_] = mutrack->normalizedChi2();
    
       }
    
    pfEvt_.nMuonpart_++;
    }

   

  // All done
  pfTree_->Fill();
  
}


bool
PFCandAnalyzer::findW(const reco::Candidate *gen){
  int nMothers = gen->numberOfMothers();                                                                 
  
  for(int iCand = 0; iCand < nMothers; ++iCand)                                                           
    {                                                                                                       
      const Candidate* currentCand = gen->mother(iCand);                                                   
      if((abs(currentCand->pdgId())==24) || abs(currentCand->pdgId())==14 ){

        if(abs(currentCand->pdgId())==14){
	  //cout<<" found mother neutrino "<<endl;
  	     return findW(currentCand);
             }

	else if(abs(currentCand->pdgId())==24){
	  // cout<<" found W "<<endl;
             
        	return true;
      }
      
      }

      else return false;
    }                                             
   return false; 
   }

bool
PFCandAnalyzer::findWbis(const reco::Candidate *gen){
  int nMothersmu = gen->numberOfMothers();                                                                 
  
  for(int iCand = 0; iCand < nMothersmu; ++iCand)                                                           
    {                                                                                                       
      const Candidate* currentCand = gen->mother(iCand);                                                   
      if((abs(currentCand->pdgId())==23) || abs(currentCand->pdgId())==13 ){

        if(abs(currentCand->pdgId())==13){
	  //cout<<" found mother neutrino "<<endl;
  	     return findWbis(currentCand);
             }

	else if(abs(currentCand->pdgId())==23){
	  // cout<<" found W "<<endl;
        	return true;
      }
      
      }

      else return false;
    }                                             
   return false; 
   }


bool
PFCandAnalyzer:: findWTris(const reco::Candidate *gen){
  int nMothersmu = gen->numberOfMothers();                                                                 
  
  for(int iCand = 0; iCand < nMothersmu; ++iCand)                                                           
    {                                                                                                       
      const Candidate* currentCand = gen->mother(iCand);                                                   
      if((abs(currentCand->pdgId())==15) || abs(currentCand->pdgId())==24 ){

        if(abs(currentCand->pdgId())==15){
	  //cout<<" found mother neutrino "<<endl;
  	     return findWbis(currentCand);
             }

	else if(abs(currentCand->pdgId())==24){
	  // cout<<" found W "<<endl;
        	return true;
      }
      
      }

      else return false;
    }                                             
   return false; 
   }


/*void HiPFCandAnalyzer::FillEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, TreePFCandEventData & tr)
{
  // General Info
  //tr.run_	  = iEvent.id().run();
  //tr.evt_	  = iEvent.id().event();
  //tr.lumi_	  = iEvent.luminosityBlock();

  if(!genOnly_&&sampleType_<10){
    // HI Event info
    edm::Handle<reco::Centrality> cent;
    iEvent.getByLabel(edm::InputTag("hiCentrality"),cent);
    Double_t hf	  = cent->EtHFhitSum();
    // Get Centrality bin
    cbins_ = getCentralityBinsFromDB(iSetup);
    tr.cent_ = cbins_->getBin(hf)*(100./cbins_->getNbins());
  }

  if (isMC_&&sampleType_<10) {
    edm::Handle<edm::GenHIEvent> mchievt;
    iEvent.getByLabel(edm::InputTag("heavyIon"),mchievt);
    tr.b_	  = mchievt->b();
    tr.npart_	  = mchievt->Npart();
    tr.ncoll_	  = mchievt->Ncoll();
  }
  }*/




void PFCandAnalyzer::beginJob()
{

  // -- trees --                                                                                                                                                                                                                        
    pfTree_ = fs->make<TTree>("pfTree","dijet tree");
    pfEvt_.SetTree(pfTree_);
    pfEvt_.doMC = doMC_;
    pfEvt_.doJets = doJets_;
    pfEvt_.SetBranches();
}

// ------------ method called once each job just after ending the event loop  ------------                                                                                                                                              
void
PFCandAnalyzer::endJob() {

  // ===== Done =====                                                                                                                                                                                                                   
  /*  if (verbosity_>=1) {
    cout << endl << "================ Ana Process Summaries =============" << endl;
    cout << " AnaJet: " << jetsrc_ << endl;
    if (refJetType_>=0) cout << " RefJet: " << refjetsrc_ << endl;
    cout << " AnaTrk: " << trksrc_ << endl;
    cout << "# HI Events : "<< numHiEvtSel_<< endl;
    cout << "# Base Events: "<< numEvtSel_ << endl;
    cout << "# Jet Events: "<< numJetEvtSel_<< endl;
  }
  */
  //beamspot->Write();
  //fhisto->Write();
  //fhisto->Close();
}

// constructors
TreePFCandEventData::TreePFCandEventData(){

}


// set branches
void TreePFCandEventData::SetBranches()
{

  
  // --event level--
   tree_->Branch("centBin",&(this->centBin),"centBin/I");
  tree_->Branch("hiNtracks",&(this->hiNtracks),"hiNtracks/I");
  tree_->Branch("hiHF",&(this->hiHF),"hiHF/F");
  tree_->Branch("hiHFplus",&(this->hiHFplus),"hiHFplus/F");
  tree_->Branch("hiHFminus",&(this->hiHFminus),"hiHFminus/F");
  tree_->Branch("hiHFplusEta4",&(this->hiHFplusEta4),"hiHFplusEta4/F");
  tree_->Branch("hiHFminusEta4",&(this->hiHFminusEta4),"hiHFminusEta4/F");
  tree_->Branch("nEvent",&(this->nEvent),"nEvent/I");
  tree_->Branch("nRun",&(this->nRun),"nRun/I");
  tree_->Branch("nLumi",&(this->nLumi),"nLumi/I");
  tree_->Branch("VertexZ",&(this->VertexZ),"VertexZ/F");

  //MET object
  tree_->Branch("recoPfMET",&(this->recoPfMET_),"recoPfMET/F");
  tree_->Branch("recoPfMETPhi",&(this->recoPfMETPhi_),"recoPfMETPhi/F");
  tree_->Branch("recoPfMETsumEt",&(this->recoPfMETsumEt_),"recoPfMETsumEt/F");
  tree_->Branch("recoPfMETSig",&(this->recoPfMETSig_),"recoPfMETSig/F");
  tree_->Branch("recoPfMETmEtSig",&(this->recoPfMETmEtSig_),"recoPfMETmEtSig/F");

  // -- particle info --
  tree_->Branch("nPFpart",&(this->nPFpart_),"nPFpart/I");
  tree_->Branch("pfId",this->pfId_,"pfId[nPFpart]/I");
  tree_->Branch("pfPt",this->pfPt_,"pfPt[nPFpart]/F");
  tree_->Branch("pfEta",this->pfEta_,"pfEta[nPFpart]/F");
  tree_->Branch("pfPhi",this->pfPhi_,"pfPhi[nPFpart]/F");
  tree_->Branch("pfTheta",this->pfTheta_,"pfTheta[nPFpart]/F");
  tree_->Branch("pfEt",this->pfEt_,"pfTheta[nPFpart]/F");
  tree_->Branch("pfCharge",this->pfCharge_,"pfCharge[nPFpart]/I");
  tree_->Branch("pfMuonpx",this->pfMuonpx_,"pfMuonpx[nPFpart]/F");
  tree_->Branch("pfMuonpy",this->pfMuonpy_,"pfMuonpy[nPFpart]/F");
  tree_->Branch("pfMuonpz",this->pfMuonpz_,"pfMuonpz[nPFpart]/F");
  tree_->Branch("pfChargepx",this->pfChargepx_,"pfChargepx[nPFpart]/F");
  tree_->Branch("pfChargepy",this->pfChargepy_,"pfChargepy[nPFpart]/F");
  tree_->Branch("pfChargepz",this->pfChargepz_,"pfChargepz[nPFpart]/F");
  //tree_->Branch("pfMuone",this->pfMuone_,"pfMuone[nPFmuon]/F");
  tree_->Branch("pfTrackdijet",this->pfTrackdijet_,"pfTrackdijet[nPFpart]/F");
  tree_->Branch("pfTrack",this->pfTrack_,"pfTrack[nPFpart]/F");
  tree_->Branch("pfTrackerMuon",this->pfTrackerMuon_,"pfTrackerMuon[nPFpart]/I");
  tree_->Branch("pfTrackhits",this->pfTrackhits_,"pfTrackhits[nPFpart]/I");
  tree_->Branch("pfdxy",this->pfdxy_,"pfdxy[nPFpart]/I");
  tree_->Branch("pfchi2",this->pfchi2_,"pfchi2[nPFpart]/I");
  tree_->Branch("pfGloMuon",this->pfGloMuon_,"pfGloMuon[nPFpart]/F");
  tree_->Branch("nGMpart",&(this->nGMpart_),"nGMpart/I");
  // -- jet info --
  if(doJets){
  tree_->Branch("njets",&(this->njets_),"njets/I");
  tree_->Branch("rawPt",this->rawPt_,"rawPt[njets]/F");
  tree_->Branch("jetPt",this->jetPt_,"jetPt[njets]/F");
  tree_->Branch("jetEta",this->jetEta_,"jetEta[njets]/F");
  tree_->Branch("jetPhi",this->jetPhi_,"jetPhi[njets]/F");
  tree_->Branch("jetMass",this->jetMass_,"jetMass[njets]/F");
  tree_->Branch("jetY",this->jetY_,"jetY[njets]/F");
  tree_->Branch("jetPU",this->jetPU_,"jetPU[njets]/F");
  tree_->Branch("jetMatchIndex",this->jetMatchIndex_,"jetMatchIndex[nPFpart]/I");  
  }

  // --tracks info--
  tree_->Branch("nTRACKpart",&(this->nTRACKpart_),"nTRACKpart/I");
  tree_->Branch("traPt",this->traPt_,"traPt[nTRACKpart]/F");
  tree_->Branch("traEta",this->traEta_,"traEta[nTRACKpart]/F");
  tree_->Branch("traPhi",this->traPhi_,"traPhi[nTRACKpart]/F");
  tree_->Branch("traQual",this->traQual_,"traQual[nTRACKpart]/I");
  tree_->Branch("traAlgo",this->traAlgo_,"traAlgo[nTRACKpart]/I");
  tree_->Branch("trahits",this->trahits_,"trahits[nTRACKpart]/I");


  //muon info
  tree_->Branch("nMuonpart",&(this->nMuonpart_),"nMuonpart/I");
  tree_->Branch("muPt",this->muPt_,"muPt[nMuonpart]/F");
  tree_->Branch("muPx",this->muPx_,"muPx[nMuonpart]/F");
  tree_->Branch("muPy",this->muPy_,"muPy[nMuonpart]/F");
  tree_->Branch("muPz",this->muPz_,"muPz[nMuonpart]/F");
  tree_->Branch("muEta",this->muEta_,"muEta[nMuonpart]/F");
  tree_->Branch("muPhi",this->muPhi_,"muPhi[nMuonpart]/F");
  tree_->Branch("muCharge",this->muCharge_,"muCharge[nMuonpart]/I");
  tree_->Branch("muTrackHits",this->muTrackHits_,"muTrackHits[nMuonpart]/I");
  tree_->Branch("muPixHits",this->muPixHits_,"muPixHits[nMuonpart]/I");
  tree_->Branch("mudxy",this->mudxy_,"mudxy[nMuonpart]/F");
  tree_->Branch("muchi2",this->muchi2_,"muchi2[nMuonpart]/I");
  tree_->Branch("chi2Track",this->chi2Track_,"chi2Track[nMuonpart]/I");
  tree_->Branch("muSeg",this->muSeg_,"muSeg[nMuonpart]/I");
  tree_->Branch("StandAloneMuon",this->StandAloneMuon_,"StandAloneMuon[nMuonpart]/I");
  tree_->Branch("stdTrack",this->stdTrack_,"stdTrack[nMuonpart]/F");
  tree_->Branch("nMuValHits",this->nMuValHits_,"nMuValHits[nMuonpart]/I");
  tree_->Branch("GMpatPt",this->GMpatPt_,"GMpatPt[nMuonpart]/F");
  tree_->Branch("mutrackPt",this->mutrackPt_,"mutrackPt[nMuonpart]/F");
  tree_->Branch("muTrackIso",this->muTrackIso_,"muTrackIso[nMuonpart]/F");
  tree_->Branch("muCaloIso",this->muCaloIso_,"muCaloIso[nMuonpart]/F");
  tree_->Branch("muEcalIso",this->muEcalIso_,"muEcalIso[nMuonpart]/F");
  tree_->Branch("muHcalIso",this->muHcalIso_,"muHcalIso[nMuonpart]/F");
  tree_->Branch("isGlobalMuon",this->isGlobalMuon_,"isGlobalMuon[nMuonpart]/O");
  tree_->Branch("isTrackerMuon",this->isTrackerMuon_,"isTrackerMuon[nMuonpart]/O");
  tree_->Branch("pass12Filter",this->pass12Filter_,"pass12Filter[nMuonpart]/O");
  tree_->Branch("pass7Filter",this->pass7Filter_,"pass7Filter[nMuonpart]/O");
  tree_->Branch("pass3Filter",this->pass3Filter_,"pass3Filter[nMuonpart]/O");
  tree_->Branch("pass12Path",this->pass12Path_,"pass12Path[nMuonpart]/O");
  tree_->Branch("pass7Path",this->pass7Path_,"pass7Path[nMuonpart]/O");
  tree_->Branch("pass3Path",this->pass3Path_,"pass3Path[nMuonpart]/O");


  

  // -- gen info --
   if(doMC){
  tree_->Branch("nGENpart",&(this->nGENpart_),"nGENpart/I");
  tree_->Branch("genStatus",this->genStatus_,"genStatus[nGENpart]/I");
  tree_->Branch("genPDGId",this->genPDGId_,"genPDGId[nGENpart]/I");
  tree_->Branch("genPt",this->genPt_,"genPt[nGENpart]/F");
  tree_->Branch("genEta",this->genEta_,"genEta[nGENpart]/F");
  tree_->Branch("genPhi",this->genPhi_,"genPhi[nGENpart]/F");
  tree_->Branch("genEt",this->genEt_,"genEt[nGENpart]/F");
  tree_->Branch("genNbMother",this->genNbMother_,"genNbMother[nGENpart]/I");
  tree_->Branch("genMotherId",this->genMotherId_,"genMotherId[nGENpart]/I");
  
  tree_->Branch("genNeuEt",&(this->genNeuEt_),"genNeuEt/F");
  tree_->Branch("genNeuEta",&(this->genNeuEta_),"genNeuEta/F");
  tree_->Branch("genNeuPhi",&(this->genNeuPhi_),"genNeuPhi/F");
  tree_->Branch("genNeuPt",&(this->genNeuPt_),"genNeuPt/F");
  
  tree_->Branch("genMuCharge",this->genMuCharge_,"genMuCharge[2]/I");
  tree_->Branch("genMuEta",this->genMuEta_,"genMuEta[2]/F");
  tree_->Branch("genMuPhi",this->genMuPhi_,"genMuPhi[2]/F");
  tree_->Branch("genMuPt",this->genMuPt_,"genMuPt[2]/F");
  tree_->Branch("genMupx",this->genMupx_,"genMupx[2]/F");
  tree_->Branch("genMupy",this->genMupy_,"genMupy[2]/F");
  tree_->Branch("nGENMu",&(this->nGENMu_),"nGENMu/I");
  tree_->Branch("nGENNeu",&(this->nGENNeu_),"nGENNeu/I");
  // tree_->Branch("nGENWout",&(this->nGENWout_),"nGENWout/I");
  /*
  tree_->Branch("genZstatus",&(this->genZstatus_),"genZstatus/I");
  tree_->Branch("genZcharge",&(this->genZcharge_),"genZcharge/I");
  tree_->Branch("genZrapidity",&(this->genZrapidity_),"genZrapidity/F");
  tree_->Branch("genZEta",&(this->genZEta_),"genZEta/F");
  tree_->Branch("genZPhi",&(this->genZPhi_),"genZPhi/F");
  tree_->Branch("genZPt",&(this->genZPt_),"genZPt/F");
  tree_->Branch("genZmass",&(this->genZmass_),"genZmass/F");
  */

  tree_->Branch("genTauPt",&(this->genTauPt_),"genTauPt/F");
  tree_->Branch("genTauEta",&(this->genTauEta_),"genTauEta/F");
  tree_->Branch("genTauPhi",&(this->genTauPhi_),"genTauPhi/F");
  tree_->Branch("genTauCharge",&(this->genTauCharge_),"genTauCharge/I");
  tree_->Branch("genTaupx",&(this->genTaupx_),"genTaupx/F");
  tree_->Branch("genTaupy",&(this->genTaupy_),"genTaupy/F");
 
         
      
  tree_->Branch("genWstatus",&(this->genWstatus_),"genWstatus/I");
  tree_->Branch("genWcharge",&(this->genWcharge_),"genWcharge/I");
  tree_->Branch("genWrapidity",&(this->genWrapidity_),"genWrapidity/F");
  tree_->Branch("genWEta",&(this->genWEta_),"genWEta/F");
  tree_->Branch("genWPhi",&(this->genWPhi_),"genWPhi/F");
  tree_->Branch("genWPt",&(this->genWPt_),"genWPt/F");
  tree_->Branch("genWmass",&(this->genWmass_),"genWmass/F");

  
  }

}
void TreePFCandEventData::Clear()
{
  // event
  genNeuEt_ =0.;
  genNeuEta_ =0;
  nPFpart_      = 0;
  njets_        = 0;
  nGENpart_     = 0;
  nGENNeu_      = 0;
  nGENMu_       = 0;
  nTRACKpart_    = 0;
  nMuonpart_    =0;
  nPFmuon_ = 0;
  nGMpart_ =0;
  nGENWout_ =0;
}




DEFINE_FWK_MODULE(PFCandAnalyzer);

