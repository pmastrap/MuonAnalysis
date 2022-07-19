import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


obj = VarParsing ('analysis')
obj.register ('inputFileList',
               '/eos/user/p/pmastrap/DYT/MuonGun/p1TeV/p1TeV.txt',
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



process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring( eos/user/p/pmastrap/DYT/MuonGun/p2TeV/p2Tev.txt
#        'file:/eos/user/p/pmastrap/DYT/MuonGun/p1TeV/step3_ext_0.root'
#    )
#)

process.source = cms.Source('PoolSource',duplicateCheckMode = cms.untracked.string("noDuplicateCheck"), fileNames = cms.untracked.vstring())

if obj.inputFileList.find(".root") < 0:
   readFileList(process.source.fileNames, obj.inputFileList,'file:')
else:
  process.source.fileNames = cms.untracked.vstring("file:"+obj.inputFileList)



process.demo = cms.EDAnalyzer('DYTthrScanTuner', psim = cms.double(1000),
                              label_track_coll = cms.InputTag("generalTracks"),
                              label_recHit_CSC  = cms.InputTag("csc2DRecHits"),
                              label_simHit_CSC  = cms.InputTag("g4SimHits","MuonCSCHits"), 
                              label_seg_CSC     = cms.InputTag("cscSegments"),
                              label_gen_coll    = cms.InputTag("genParticles"),
                              label_muon        = cms.InputTag("muons")
                              #open = cms.string('recreate'),
                              #out = cms.string('alessio.root')

)

process.TFileService = cms.Service("TFileService", fileName = cms.string('alessio.root'))
process.p = cms.Path(process.demo)
