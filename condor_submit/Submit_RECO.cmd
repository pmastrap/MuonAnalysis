Universe                = vanilla
Environment             = CONDORJOBID=$(Process)
when_to_transfer_output = ON_EXIT
transfer_output_files   = ""
transfer_input_files    = /afs/cern.ch/user/p/pmastrap/private/CSC/CMSSW_12_2_1/src/MuonGun/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py 
should_transfer_files   = YES
executable              = run_the_matrix_RECO.sh
arguments               = $(ProcId)
output                  = logs/$Fn(executable).$(ClusterId).$(ProcId).out
error                   = logs/$Fn(executable).$(ClusterId).$(ProcId).err
log                     = logs/$Fn(executable).$(ClusterId).$(ProcId).log
+JobFlavour             = "espresso"
queue 2 
