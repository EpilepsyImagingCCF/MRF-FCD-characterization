MNI_template=/path/to/MNI152_T1_1mm_brain.nii

target_dir=/path/to/target/

transform_type=s

# register the MRF T1w brain file into the MNI space
antsRegistrationSyN.sh -d 3 \
	-m ${target_dir}/T1w_data_brain.nii \
        -f ${MNI_template} \
	-t ${transform_type} \
	-n 3 \
        -o ${target_dir}/Reg_MRF_2_MNI_${transform_type}_

# Apply the transformation matrix to the MRF T1w, T1 and T2 brain image
antsApplyTransforms -d 3 \
        -i ${target_dir}/T1w_data_brain.nii \
        -o ${target_dir}/n_syN_T1w_data_brain_Warped.nii \
        -r ${MNI_template} \
        -t ${target_dir}/Reg_MRF_2_MNI_${transform_type}_1Warp.nii.gz \
        -t ${target_dir}/Reg_MRF_2_MNI_${transform_type}_0GenericAffine.mat

antsApplyTransforms -d 3 \
        -i ${target_dir}/T1_data_brain.nii \
        -o ${target_dir}/n_syN_T1_data_brain_Warped.nii \
        -r ${MNI_template} \
        -t ${target_dir}/Reg_MRF_2_MNI_${transform_type}_1Warp.nii.gz \
        -t ${target_dir}/Reg_MRF_2_MNI_${transform_type}_0GenericAffine.mat

antsApplyTransforms -d 3 \
        -i ${target_dir}/T2_data_brain.nii \
        -o ${target_dir}/n_syN_T2_data_brain_Warped.nii \
        -r ${MNI_template} \
        -t ${target_dir}/Reg_MRF_2_MNI_${transform_type}_1Warp.nii.gz \
        -t ${target_dir}/Reg_MRF_2_MNI_${transform_type}_0GenericAffine.mat

# Apply the transformation matrix to the ROI
antsApplyTransforms -d 3 \
        -i ${target_dir}/FCD_ROI.nii \
        -o ${target_dir}/n_FCD_ROI.nii \
        -r ${MNI_template} \
        -t ${target_dir}/Reg_MRF_2_MNI_${transform_type}_1Warp.nii.gz \
        -t ${target_dir}/Reg_MRF_2_MNI_${transform_type}_0GenericAffine.mat

fslmaths ${target_dir}/n_FCD_ROI.nii -thr 0.5 -bin ${target_dir}/n_FCD_ROI_bin.nii
