import FWCore.ParameterSet.Config as cms

l1JetRecoTree = cms.EDAnalyzer("L1JetRecoTreeProducer",
  pfJetToken              = cms.untracked.InputTag("ak4PFJetsCHS"),
  jecToken                = cms.untracked.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector"),
  caloJetToken            = cms.untracked.InputTag("ak4CaloJets"),
  caloJetIDToken          = cms.untracked.InputTag("ak4JetID"),
  maxJet                  = cms.uint32(20),
  jetptThreshold          = cms.double(30),
  jetetaMax               = cms.double(2.5),
  pfMetToken              = cms.untracked.InputTag("pfMetT1"),
  caloMetToken            = cms.untracked.InputTag("caloMet"),
  caloMetBEToken          = cms.untracked.InputTag("caloMetBE")
)


