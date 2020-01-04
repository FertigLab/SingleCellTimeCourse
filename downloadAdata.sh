#!/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>log.out 2>&1
## get folder names and store in array
files=( $(aws s3 ls s3://fertig-lab-bucket/projects/Short_Time_Course_Single_Cell/ | grep output | awk '{print $2}') )
#folders=( ${folders[@]:2} )
echo ${files[3]}

# download fileGroup files
for file in ${files[@]}; do
	#download one file SXX_L00X.txt
	echo downloading file $file
	fd=$(echo $file | tr -d '/' | cat)
	aws s3 cp s3://fertig-lab-bucket/projects/Short_Time_Course_Single_Cell/${file}counts_unfiltered/adata.h5ad ./adataObjs/${fd}_adata.h5ad
done
