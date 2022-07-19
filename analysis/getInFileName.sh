#! /bin/sh
## MAKE A TXT FILE WITH ALL THE INPUT NTUPLES

INPUT_FOLDER_I=/eos/user/p/pmastrap/DYT/MuonGun/p2TeV/
#INPUT_FOLDER_II=/eos/user/p/pmastrap/DYT/MuonGun/p02TeV/

cd ${INPUT_FOLDER_I}

for d in */ ; do

    echo ${d///} #"$d" 
    ls -ltr | awk '/step3/ {print "/eos/user/p/pmastrap/DYT/MuonGun/p2TeV/" $9;}' > p2TeV.txt
    
done

#cd ${INPUT_FOLDER_II}

#for d in */ ; do
    
#    echo ${d///} #"$d" 
#    ls -ltr | awk '/step3/ {print "/eos/user/p/pmastrap/DYT/MuonGun/p02TeV/" $9;}' > p02TeV.txt
    
#done


