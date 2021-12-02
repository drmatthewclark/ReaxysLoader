#!/bin/bash
# download latest rmc dataset

echo -n "reaxys "
source update.sh
source ./ReaxysLoader/credentials.py

update $reaxys_loc $reaxys_name

cd ${latest}
eval "$(conda shell.bash hook)"
conda activate standard
time python ../ReaxysLoader/load_all.py
../../fix-perms

cd ..

del() {
  shift
  for d in $*; do
        echo 'removing dataset' $d
        rm -r "${d}"
  done

}

dirs=`ls -dc 2[0-9]*`
del ${dirs}
