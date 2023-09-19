#! /bin/bash

# This script should be run AFTER Final_Recon_All.sh

N=5
for dir in /mnt/udata/shared/database/meg-ieeg-UNMC/derivatives/*/ses-*op; do
	(
		export SUBJECTS_DIR=${dir}

		# check for sample folder
		if [ -d ${SUBJECTS_DIR}/sample ]; then

			# check for bem folder
			if [ -d ${SUBJECTS_DIR}/sample/bem ]; then
				echo "${SUBJECTS_DIR} already processed."
			else
				echo "Starting task ${SUBJECTS_DIR}."
				SUBJECT=sample

				# Create BEM surfaces using watershed algorithm
				mne watershed_bem -s $SUBJECT -d $SUBJECTS_DIR -o -b ws.mgz

				# Create high resolution head surfaces for coordinate alignment
				mne make_scalp_surfaces -s $SUBJECT -d $SUBJECTS_DIR -o -f

			fi
		else
			echo "${SUBJECTS_DIR}/sample/ does not exist."
		fi
	) &

	# Execute N jobs in parallel
	if [[ 'jobs -r -p | wc -l' -ge $N ]]; then
		wait -n
	fi

done
wait
echo "All Done!"
