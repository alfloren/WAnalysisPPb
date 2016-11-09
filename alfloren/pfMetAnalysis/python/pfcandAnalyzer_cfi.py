import FWCore.ParameterSet.Config as cms

pfcandAnalyzer = cms.EDAnalyzer('PFCandAnalyzer',
                                pfCandidateLabel = cms.InputTag("particleFlowTmp"),
                                jetLabel = cms.InputTag("ak5PFpatJets"),
                                genLabel = cms.InputTag("genParticles"),
                                recoPfMETLabel = cms.InputTag("pfMet"),
                                VtxLabel = cms.InputTag("offlinePrimaryVertices"),
                                # debug
                                pfPtMin = cms.double(0.5),
                                genPtMin = cms.double(0.5),
                                jetPtMin = cms.double(20.0),                                
                                verbosity = cms.untracked.int32(0),
                                doJets = cms.untracked.bool(False),
                                doMC = cms.untracked.bool(True),
                                skipCharged = cms.untracked.bool(False)
                                )

