#! /bin/bash

# After this function, run Final_Create_Surfaces.sh

export FREESURFER_HOME=/usr/local/freesurfer/7.2.0
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

N=5

inDir=/mnt/udata/shared/database/meg-ieeg-UNMC/rawdata
outDir=/mnt/udata/shared/database/meg-ieeg-UNMC/derivatives

for dir in $inDir/*/ses-*; do
	(
		export SUBJECTS_DIR=${dir}
		#echo "${SUBJECTS_DIR}"

		sub=$(echo $dir| cut -d'/' -f 8)
		op=$(echo $dir| cut -d'/' -f 9)

		outDir_full=${outDir}/${sub}/${op}
		#echo "${outDir_full}"

		# Check for anatomy folder
		if [ -d "${SUBJECTS_DIR}/anat" ]; then
		
			# Check for sample folder or test.txt
			if [ -d "${SUBJECTS_DIR}/sample/" ] || [ -d "${outDir_full}/sample/" ] || [ -d "${SUBJECTS_DIR}/test.txt" ]; then
				echo "${SUBJECTS_DIR}/sample/ exists."
			else
				touch ${SUBJECTS_DIR}/test.txt
				my_subject=sample
				my_NIfTI=${SUBJECTS_DIR}/anat/*T1w.nii.gz
				recon-all -i $my_NIfTI -s $my_subject -all
				rm ${SUBJECTS_DIR}/test.txt
			fi
		else
			echo "${SUBJECTS_DIR}/anat/ does not exist."
		fi

		# Move sample folder to derivatives folder
		if [ -d "${SUBJECTS_DIR}/sample/" ]; then
			if [ -d "${outDir_full}" ]; then
				mv ${SUBJECTS_DIR}/sample ${outDir_full}
			else
				mkdir -p ${outDir_full}
				mv ${SUBJECTS_DIR}/sample ${outDir_full}
			fi
		fi
	)&

	# Execute N jobs in parallel
	if [[ `jobs -r -p | wc -l` -ge $N ]]; then
		wait -n
	fi
done
wait
echo "All Done!"
