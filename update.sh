#!/bin/bash 

# function to download update 
update() {

data=$1     # example: rmc-ff-download
profile=$2  # example: rmc

list=`aws  --profile ${profile} s3 ls ${data}/ | sed 's/[^0-9]*//g' | sort -n -r | head -n 2`

for release in $list; do
	if [ -d ${release} ]; then
	       echo "${release} already downloaded"
	       exit 0 # already downloaded
	fi

	aws --profile ${profile}  s3 cp s3://${data}/${release} ./${release} --recursive --no-progress
	ecode=$?

	if [ $ecode = 0 ]; then 
		echo "...${release} downloaded"
		return 0 # exit if success
	else
		echo "error downloading ${release}"
		rm -r "${release}"
	fi
done
}
