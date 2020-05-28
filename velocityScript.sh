#!/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>log.out 2>&1
## get folder names and store in array
files=( $(aws s3 ls s3://fertig-lab-bucket/projects/Short_Time_Course_Single_Cell/fileGroups/ | awk '{print $4}') )
#folders=( ${folders[@]:2} )
#echo ${folders[@]}

# download fileGroup files
for file in ${files[@]}; do
	#download one file SXX_L00X.txt
	echo downloading file $file
	aws s3 cp s3://fertig-lab-bucket/projects/Short_Time_Course_Single_Cell/fileGroups/${file} ./
	fastqFiles=( $(awk '{print $1}' $file))
	folderName=$(echo $file | tr -d '.txt' | cat)
	echo making folder $folderName
	mkdir $folderName
	rfq=()
	for fq in ${fastqFiles[@]}; do
		fq=$(echo $fq | tr -d '\040\011\012\015' | cat)
		aws s3 cp s3://fertig-lab-bucket/projects/Short_Time_Course_Single_Cell/${fq} ${folderName}/
		if echo $fq | grep _R; then
			rfq+=( ${folderName}/${fq} )
		fi
    	done
	echo running kallisto bustools for $folderName
	kb count -i homo_sapiens/index.idx -o output_$folderName -x 10xv2 -t 36 -g homo_sapiens/transcripts_to_genes.txt -c1 homo_sapiens/cdna_transcripts_to_capture.txt -c2 homo_sapiens/intron_transcripts_to_capture.txt --lamanno --h5ad ${rfq[@]}
	echo removing folder $folderName
	rm -r $folderName
	echo removing file $file
	rm $file
	echo uploading the output
	aws s3 cp ./output_${folderName}/ s3://fertig-lab-bucket/projects/Short_Time_Course_Single_Cell/output_${folderName}/ --recursive
	echo removing output
	rm -r output_${folderName}
done
