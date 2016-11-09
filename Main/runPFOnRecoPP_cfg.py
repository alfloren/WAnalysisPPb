### -------------------------------------------------------------------
### VarParsing allows one to specify certain parameters in the command line
### e.g.
### cmsRun testElectronSequence_cfg.py print maxEvents=10
### -------------------------------------------------------------------

import FWCore.ParameterSet.Config as cms
#import FWCore.ParameterSet.VarParsing as VarParsing
#import os 

process = cms.Process("PFRECO")


process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load("CommonTools.UtilAlgos.TFileService_cfi")

#global tags for conditions data: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_P_V43D::All'
##################################################################################

# setup 'standard'  options
#options = VarParsing.VarParsing ('standard')

# setup any defaults you want
#options.output = 'test_out.root'
#options.files = [
    #'/store/relval/CMSSW_4_3_0/RelValHydjetQ_MinBias_2760GeV/GEN-SIM-RECO/STARTHI43_V1-v1/0125/8C6E49D1-11A1-E011-B2A6-003048678F8E.root'
    #'/store/relval/CMSSW_3_9_0/RelValPyquen_ZeemumuJets_pt10_2760GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_39Y_V2-v1/0052/4AFA1DFA-1AD9-DF11-BBEB-003048679046.root'
 #   ] 
#options.maxEvents = 1 

# get and parse the command line arguments
#options.parseArguments()


##################################################################################
# Some Services
	   
#process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
#                                        ignoreTotal=cms.untracked.int32(0),
#                                        oncePerEventMode = cms.untracked.bool(False)
#                                        )
#
#process.Timing = cms.Service("Timing")

##################################################################################
# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

# pile up rejection
process.load('Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cff')

# Centrality for pPb
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')

#process.load('HeavyIonsAnalysis.Configuration.CommonFunctions_cff')
#overrideCentrality(process)
process.GlobalTag.toGet = cms.VPSet(
   cms.PSet(record = cms.string("HeavyIonRcd"),
                   tag = cms.string("CentralityTable_HFtrunc100_PA2012B_v538x02_offline"),
                   connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS"),
                   label = cms.untracked.string("HFtowersTrunc")
                   )
    )



process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowersTrunc"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("pACentrality"),
    pPbRunFlip = cms.untracked.uint32(211313)
    )


from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cff import *
onia2MuMuPAT(process, GlobalTag=process.GlobalTag.globaltag, MC=False, HLT="HLT", Filter=False)

process.patMuonSequence = cms.Sequence(
    process.PAcollisionEventSelection *
    process.pileupVertexFilterCutGplus * 
    process.pACentrality_step *
    process.patMuonsWithTriggerSequence
    )

# W part
process.goodPatMuons = cms.EDFilter("PATMuonSelector",
                                    src = cms.InputTag("patMuonsWithTrigger"),
                                    cut = cms.string("pt>15."),
                                    #cut = cms.string(""),
                                    filter = cms.bool(True)
                                    )

# Input Source
process.source.fileNames = cms.untracked.vstring(
    '/store/hidata/HIRun2013/PAMuon/RECO/PromptReco-v1/000/210/498/00000/C08874D8-F364-E211-9A31-BCAEC518FF67.root'
)

#process.source.skipEvents = cms.untracked.uint32(600)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(
    -1)
)
    
process.load("alfloren/pfMetAnalysis/pfcandAnalyzerData_cfi")
process.pfcandAnalyzerData.pfCandidateLabel = "particleFlow"
process.pfcandAnalyzerData.genLabel = "GenParticles"
#process.ak5PFJetAnalyzer.

process.load("MitHig.PixelTrackletAnalyzer.METAnalyzer_cff")
# makes the PF Towers
process.load("RecoHI.HiJetAlgos.ParticleTowerProducer_cff")
process.particleTowerProducer.src = "particleFlow"
# puts them in an ntuple
process.load("alfloren.pfMetAnalysis.pftowerAnalyzer_cff")

process.anaMET.METSrc = 'genMetTrue'

process.hltSglMu = cms.EDFilter("HLTHighLevel",
                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                 HLTPaths = cms.vstring("HLT_PAMu12_v*"),
                 eventSetupPathsKey = cms.string(''),
                 andOr = cms.bool(True),
                 throw = cms.bool(False)
)

process.load("CmsHi.JetAnalysis.PatAna_cff")
# Careful, you need a different version for MC!!!
#process.load("CmsHi.JetAnalysis.PatAna_MC_cff")


# change to proper JEC
from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_pp2760
overrideJEC_pp2760(process)
process.ak5PFcorr.payload = cms.string('AK5PF_generalTracks')

process.load("RecoHI.HiJetAlgos.HiRecoPFJets_cff")
process.akPu3PFJets = process.ak3PFJets.clone(
    rParam = cms.double(0.3),
    src = cms.InputTag("particleTowerProducer"),
    jetType = cms.string('BasicJet')
)

process.patJetPath = cms.Path(process.ak5PFpatSequence*
                              process.particleTowerProducer*
                              process.akPu3PFJets*
                              process.akPu3PFpatSequence
                              )


#AK5PF_generalTracks
process.load("CmsHi.JetAnalysis.inclusiveJetAnalyzer_cff")


process.ak5PFJetAnalyzer.matchTag = cms.untracked.InputTag("akPu3PFpatJets")
process.ak5PFJetAnalyzer.trackTag = cms.InputTag("generalTracks")
# these are default, but reset them for MC
#process.ak5PFJetAnalyzer.isMC = cms.untracked.bool(False)
#process.ak5PFJetAnalyzer.fillGenJets = cms.untracked.bool(False)

process.akPu3PFJetAnalyzer = process.ak5PFJetAnalyzer.clone(
    jetTag = cms.InputTag("akPu3PFpatJets"),
    matchTag = cms.untracked.InputTag("ak5PFpatJets")
    )


process.ntuples = cms.Path(
    process.patMuonSequence*
    #process.hltSglMu *
    process.goodPatMuons*
    process.pfcandAnalyzerData*
    process.ak5PFJetAnalyzer*
    process.akPu3PFJetAnalyzer
    #process.pftowerAna
    #process.anaMET
    )

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("WppbRunsAllrunminus7NewCentrality.root"))



#process.out_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(process.patJetPath,process.ntuples)
