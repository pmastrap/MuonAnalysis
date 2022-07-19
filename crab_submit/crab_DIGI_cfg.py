from CRABClient.UserUtilities import config
config = config()

config.General.transferOutputs = True
config.General.transferLogs = False 

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = './MuonGun/step2_DIGI_L1_DIGI2RAW_HLT.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/user/pmastrap/MuonGun_DYT' #% (getUsernameFromSiteDB())
config.Data.publication = True

config.Site.storageSite = 'T2_BE_UCL'
#config.Site.whitelist = 'T2_BE_UCL'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from http.client import HTTPException 
    
    samples = {'SingleMuonP500_DIGI_l1_DIG2RAW_HLT':  '/SingleMuonP500_pythia8/pmastrap-GEN-SIM_CMSSW_12_4_0_phase1_2018_realistic_DYTStudies-669d61c53cd506b0d351a5c7089fbeda/USER',
               'SingleMuonP1000_DIGI_l1_DIG2RAW_HLT': '/SingleMuonP1000_pythia8/pmastrap-GEN-SIM_CMSSW_12_4_0_phase1_2018_realistic_DYTStudies-881c117e1b5571ec9b6059d8af63cf5b/USER',
	       'SingleMuonP2000_DIGI_l1_DIG2RAW_HLT': '/SingleMuonP2000_pythia8/pmastrap-GEN-SIM_CMSSW_12_4_0_phase1_2018_realistic_DYTStudies-9feef9d1876f805c6687ae4020930f59/USER',
               'SingleMuonP3000_DIGI_l1_DIG2RAW_HLT': '/SingleMuonP3000_pythia8/pmastrap-GEN-SIM_CMSSW_12_4_0_phase1_2018_realistic_DYTStudies-8bb984c2de5385c4ba8a6f8ea8f5b632/USER',
               'SingleMuonP4000_DIGI_l1_DIG2RAW_HLT': '/SingleMuonP4000_pythia8/pmastrap-GEN-SIM_CMSSW_12_4_0_phase1_2018_realistic_DYTStudies-c0f5fa9b6b3c93231ebc068f9c856527/USER',
               'SingleMuonP5000_DIGI_l1_DIG2RAW_HLT': '/SingleMuonP5000_pythia8/pmastrap-GEN-SIM_CMSSW_12_4_0_phase1_2018_realistic_DYTStudies-c62a98fada06d84aefdbb21736eca7f0/USER'}    
              
    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects_DYT_MuonPGun_DIGI'

    def submit(config):
        try:
            crabCommand('submit', config = config, dryrun=False)
        except HTTPException as hte:
            print("Failed submitting task: ", (hte.headers))
        except ClientException as cle:
            print("Failed submitting task: ", (cle))

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
    for sample, gen_dataset in samples.items():
         config.General.requestName = sample
         config.Data.inputDataset = gen_dataset
         config.Data.inputDBS = 'phys03'
         #config.Data.outputPrimaryDataset = sample
         config.Data.outputDatasetTag = 'DIGI_CMSSW_10_6_X_phase1_2018_realistic_DYTStudies'
         submit(config)
