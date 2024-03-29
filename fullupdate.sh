#!/bin/bash
# download latest rmc dataset

loader="ReaxysLoader"

echo "reaxys "
source update.sh
source ./${loader}/credentials.py

update ${download} ${dataset}

if [ "${success}" = "no" ]; then
	echo "loading ended"
	exit 1
fi

cd ${release}
pwd
eval "$(conda shell.bash hook)"
conda activate standard
time python ../${loader}/load_all.py
../../fix-perms

cd ..

del() {
  shift
  for d in $*; do
	if [ ! "${d}" == "${release}" ]; then
          echo 'removing dataset' $d
          rm -r "${d}"
	fi
  done

}

dirs=`ls -dc 2[0-9]*`
del ${dirs}
