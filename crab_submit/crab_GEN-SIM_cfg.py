import CRABClient
from CRABClient.UserUtilities import config
from multiprocessing import Process

config = config()

config.General.transferOutputs = True
config.General.transferLogs = False 

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = './MuonGun/SingleMuPvar_EtaLessThan4_cfi_GEN_SIM.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = "EventBased"
config.Data.totalUnits = 2000000
config.Data.unitsPerJob = 1000

config.Data.outLFNDirBase = '/store/user/pmastrap/MuonGun_DYT_prova' #% (getUsernameFromSiteDB())
config.Data.publication = True

config.Site.storageSite = 'T2_BE_UCL'

if __name__ == '__main__':

	from CRABAPI.RawCommand import crabCommand
	from CRABClient.ClientExceptions import ClientException
	from http.client import HTTPException

	samples = {'SingleMuonP500_pythia8': 500,
                   'SingleMuonP1000_pythia8': 1000,
		   'SingleMuonP2000_pythia8': 2000,
		   'SingleMuonP3000_pythia8': 3000,
                   'SingleMuonP4000_pythia8': 4000,
                   'SingleMuonP5000_pythia8': 5000
                   }
	# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
	# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
	config.General.workArea = 'crab_projects_DYT_MuonPGun_prova'

	def submit(config):
		try:
		    crabCommand('submit', config = config)
		except HTTPException as hte:
		    print("Failed submitting task: ", (hte.headers))
		except ClientException as cle:
		    print("Failed submitting task: ", (cle))

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################
for sample, p in samples.items():
         config.General.requestName = sample
         config.JobType.pyCfgParams = ['p={0}'.format(p)]
         config.Data.outputPrimaryDataset = sample
         config.Data.outputDatasetTag = 'GEN-SIM_CMSSW_12_4_0_phase1_2018_realistic_prova'
         #submit(config)
         p = Process(target=submit, args=(config,))
         p.start()
         p.join()
