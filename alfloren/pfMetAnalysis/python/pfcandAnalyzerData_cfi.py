import FWCore.ParameterSet.Config as cms

pfcandAnalyzerData = cms.EDAnalyzer('PFCandAnalyzerData',
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
                                doJets = cms.untracked.bool(True),
                                skipCharged = cms.untracked.bool(False)
                                )

