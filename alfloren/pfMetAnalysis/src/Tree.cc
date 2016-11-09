#include "alfloren/pfMetAnalysis/interface/Tree.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"
using namespace std;


// constructors
Tree::Tree(){

}


// set branches
void Tree::SetBranches()
{
  // --event level--

  // -- particle info --
  tree_->Branch("nPFpart",&(this->nPFpart_),"nPFpart/I");
  tree_->Branch("pfId",this->pfId_,"pfId[nPFpart]/I");
  tree_->Branch("pfPt",this->pfPt_,"pfPt[nPFpart]/F");
  tree_->Branch("pfEta",this->pfEta_,"pfEta[nPFpart]/F");
  tree_->Branch("pfPhi",this->pfPhi_,"pfPhi[nPFpart]/F");

  // -- jet info --
  tree_->Branch("njets",&(this->njets_),"njets/I");
  tree_->Branch("jetPt",this->jetPt_,"jetPt[njets]/F");
  tree_->Branch("jetEta",this->jetEta_,"jetEta[njets]/F");
  tree_->Branch("jetPhi",this->jetPhi_,"jetPhi[njets]/F");

  // -- gen info --
  tree_->Branch("nGENpart",&(this->nGENpart_),"nGENpart/I");
  tree_->Branch("genPDGId",this->genPDGId_,"genPDGId[nGENpart]/I");
  tree_->Branch("genPt",this->genPt_,"genPt[nGENpart]/F");
  tree_->Branch("genEta",this->genEta_,"genEta[nGENpart]/F");
  tree_->Branch("genPhi",this->genPhi_,"genPhi[nGENpart]/F");


}
void Tree::Clear()
{
  // event

  nPFpart_      = 0;
  njets_        = 0;
  nGENpart_     = 0;
}

