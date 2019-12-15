#!/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>log.out 2>&1
## get folder names and store in array
folders=( $(aws s3 ls s3://fertig-lab-bucket/projects/Evanthia_MDSC_Analysis/ | awk '{print $2}' | grep ET01) )
folders=( ${folders[@]:2} )
#echo ${folders[@]}

# start loop to perform kallisto bus
for name in ${folders[@]}; do
	echo creating directory $name
	mkdir $name
	# list files on aws
	files=( $(aws s3 ls s3://fertig-lab-bucket/projects/Evanthia_MDSC_Analysis/$name | awk '{print $4}' | grep _R) )
	for file in ${files[@]}; do
		#echo $name$file
		#download files from aws
		aws s3 cp s3://fertig-lab-bucket/projects/Evanthia_MDSC_Analysis/$name$file $name
	done
	echo Running kallisto bustools
	kb count -i mouse_velocity_index2/index.idx -o output_$name -x 10xv2 -t 36 -g mouse_velocity_index2/transcripts_to_genes.txt \
	 -c1 mouse_velocity_index2/cdna_transcripts_to_capture.txt -c2 mouse_velocity_index2/intron_transcripts_to_capture.txt --lamanno \
	  --h5ad $name${files[0]} $name${files[1]} $name${files[2]} $name${files[3]}
	echo Generated spliced and unspliced data for $name
	echo deleting directory $name
	rm -r $name
	echo uploading output folder to aws
	aws s3 cp output_$name s3://fertig-lab-bucket/projects/Evanthia_MDSC_Analysis/output_$name --recursive
	echo deleting output folder
	rm -r output_$name
done
