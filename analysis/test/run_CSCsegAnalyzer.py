##########################################
##Author: Paola Mastrapasqua UCLouvain(BE)
##########################################

import FWCore.ParameterSet.Config as cms
import subprocess
import os

from FWCore.ParameterSet.VarParsing import VarParsing

obj = VarParsing ('analysis')
obj.register ('inputFileList',
               #os.environ.get('CMSSW_BASE')+'/src/MuonGun/Results/input_list.txt',
               '/eos/user/p/pmastrap/DYT/MuonGun/p02TeV/p02TeV.txt',
               VarParsing.multiplicity.singleton,
               VarParsing.varType.string,
               info="input files list")
# get and parse the command line arguments
obj.parseArguments()

process = cms.Process("CMSPOStiming")

def readFileList(fileList, inputFileName, fileNamePrefix):
    inputFile = open(inputFileName, 'r')
    for name in inputFile:
        fileList.extend([ fileNamePrefix + name ])
        print('Added:'+ str(fileNamePrefix) + str(name))
    inputFile.close()


process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')


#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag.globaltag = "94X_dataRun2_v11"

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.source = cms.Source('PoolSource',duplicateCheckMode = cms.untracked.string("noDuplicateCheck"), fileNames = cms.untracked.vstring())

if obj.inputFileList.find(".root") < 0:
   readFileList(process.source.fileNames, obj.inputFileList,'file:')
else:
  process.source.fileNames = cms.untracked.vstring("file:"+obj.inputFileList)

#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#    '/store/data/Run2017C/SingleMuon/AOD/17Nov2017-v1/70000/E2E0CA7E-7ADA-E711-9AC5-02163E019D1B.root'
#    )
#)

from RecoMuon.TrackingTools.MuonSegmentMatcher_cff import *

process.segmentNtupleFiller = cms.EDAnalyzer("CSCsegAnalyzer",
    MuonSegmentMatcher,

    Muons = cms.untracked.InputTag("muons"),
    TKtracks = cms.untracked.InputTag("generalTracks"),
    PrimaryVertex = cms.untracked.InputTag("offlinePrimaryVertices"),
    CSCg4SimHits = cms.untracked.InputTag("g4SimHits","MuonCSCHits"),
    g4genParticles = cms.untracked.InputTag("genParticles"),    

    open = cms.string('recreate'),
    out = cms.string('prova.root'),
)

process.jobPath = cms.Path(process.segmentNtupleFiller)

