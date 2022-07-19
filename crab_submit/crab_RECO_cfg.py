from CRABClient.UserUtilities import config
config = config()

config.General.transferOutputs = True
config.General.transferLogs = False 

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = './MuonGun/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/user/pmastrap/MuonGun_CSC_ME11_NEW/' 
#config.Data.publication = True

config.Site.storageSite = 'T2_BE_UCL'
# UCL added ot black list because doesn't have IB
# config.Site.blacklist   = ['T2_BE_UCL']
# config.Site.whitelist   = ['T2_FR_*','T2_CH_CERN','T2_IT_*']

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from http.client import HTTPException
    
    samples = {#'SingleMuonP500_RECO' : 
               'SingleMuonP1000_RECO': '/SingleMuonP1000_pythia8/pmastrap-DIGI_CMSSW_10_6_X_phase1_2018_realistic_CSCME11Studies-6f6c9f7a8c605514dd44d2880da1f919/USER',
               'SingleMuonP2000_RECO': '/SingleMuonP2000_pythia8/pmastrap-DIGI_CMSSW_10_6_X_phase1_2018_realistic_CSCME11Studies-6f6c9f7a8c605514dd44d2880da1f919/USER'  
               #'SingleMuonP3000_RECO':
               #'SingleMuonP4000_RECO':
               #'SingleMuonP5000_RECO':
               }
    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'ppp'

    def submit(config):
        try:
            crabCommand('submit', config = config, dryrun=True)
        except HTTPException as hte:
            print("Failed submitting task: ", (hte.headers))
        except ClientException as cle:
            print("Failed submitting task: ",  (cle))

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
    for sample, gen_dataset in samples.items():
         config.General.requestName = sample#+'_APE'
         #config.JobType.pyCfgParams = ['runAPEs=True']
         config.Data.inputDataset = gen_dataset
         config.Data.inputDBS = 'phys03'
         #config.Data.outputPrimaryDataset = sample
         config.Data.outputDatasetTag = 'ppp'
         submit(config)

    #for sample, gen_dataset in samples.items():
    #     config.General.requestName = sample+'_noAPE'
    #     config.JobType.pyCfgParams = ['runAPEs=False']
    #     config.Data.inputDataset = gen_dataset
    #     config.Data.inputDBS = 'phys03'
         #config.Data.outputPrimaryDataset = sample
         #config.Data.outputDatasetTag = 'CMSSW_10_6_X_phase1_2018_realistic_step3_April2019'
    #     submit(config)
