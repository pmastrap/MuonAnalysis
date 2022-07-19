#! /bin/sh
####################################
#        LaunchOnFarm Script       #
#     Loic.quertenmont@cern.ch     #
#            April 2010            #
####################################

ProcId=$1
echo ${ProcId}
# Script stop when a line fail
#set -e
#source /cvmfs/cms.cern.ch/cmsset_default.sh
LAUNCH_FOLDER=/afs/cern.ch/user/p/pmastrap/private/CSC/CMSSW_12_2_1/src/MuonGun
export SCRAM_ARCH=slc7_amd64_gcc900

cd $LAUNCH_FOLDER

echo $PWD 
 
if [[ ! -d 'Results' ]]; then
	mkdir Results
fi

eval `scramv1 runtime -sh`

cd ..

#echo $PWD
#ls -lrt

#cp ${LAUNCH_FOLDER}/Results/step2_${ProcId}.root .
cp  /eos/user/p/pmastrap/DYT/MuonGun/p2TeV/step2_${ProcId}.root .

ls ${LAUNCH_FOLDER}
ls -lrt
 

printf 'Running RECO step old DyT....\n'
cmsRun ${LAUNCH_FOLDER}/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py inputFile=step2_${ProcId}.root 
printf '....RECO old DyT finished\n'


mv step3.root step3_ext_${ProcId}.root 
mv step3_ext*.root /eos/user/p/pmastrap/DYT/MuonGun/p2TeV 
rm step2_${ProcId}.root
