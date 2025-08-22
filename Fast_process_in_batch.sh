folder=$1

for sub in ${folder}/*
do

filename=`basename ${sub}/n_syN_T1w_data_brain_Warped.nii`

echo "GM/WM/CSF seg"
fast -g ${sub}/${filename}.nii

done