//
// Class to represent event tree data
//

#ifndef PFCANDEVENTDATA_H
#define PFCANDEVENTDATA_H
#include <vector>
#include "TTree.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <vector>

const Int_t MAXPARTICLE = 100000;
//
// DiJet ana Event Data Tree definition
//
class Tree
{
  public:
    // ===== Class Methods =====
    void SetDefaults();
    Tree();
    void SetTree(TTree * t) { tree_=t; }
    void SetBranches();
    void Clear();


    Float_t		    jdphi_;


    // -- particle info --
    Int_t		    nPFpart_, nGENpart_, njets_;
    Int_t	  	    pfId_[MAXPARTICLE], genPDGId_[MAXPARTICLE];
    Float_t	  	    pfPt_[MAXPARTICLE], genPt_[MAXPARTICLE],  jetPt_[MAXPARTICLE];
    Float_t	  	    pfEta_[MAXPARTICLE], genEta_[MAXPARTICLE],  jetEta_[MAXPARTICLE];
    Float_t	  	    pfPhi_[MAXPARTICLE], genPhi_[MAXPARTICLE],  jetPhi_[MAXPARTICLE];


  private:
    TTree*                 tree_;
};
#endif // PFCANDEVENTDATA_H
