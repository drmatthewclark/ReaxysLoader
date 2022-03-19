#!/bin/bash 

# function to download update 
update() {

data=$1
profile=$2

list=`aws  --profile ${profile} s3 ls ${data}/ | sed 's/[^0-9]*//g' | sort -n -r | head -n 2`
latest=`ls -dc 2[0-9]* | tail -n 1`
echo "latest is " $latest

success="no"
for release in $list; do
	
	if [ -d ${release} ]; then
	       echo "${release} already downloaded"

	elif [ "$release" -gt "$latest" ]; then	

	  aws --profile ${profile}  s3 cp s3://${data}/${release} ./${release} --recursive --no-progress
	  ecode=$?

	  if [ $ecode = 0 ]; then 
		success="yes"
		echo "...${release} downloaded"
		return 0 # exit if success
 	  else
		echo "error downloading ${release}"
		rm -r "${release}"
	  fi
	fi
done
}
