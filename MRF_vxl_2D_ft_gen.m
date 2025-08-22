clear all; close all;

toolbox_path = '/path/ to/ toolbox';

addpath(fullfile(toolbox_path, 'NIfTI_20140122'))

pts_folder = '/path/ to/ patient folder';
controls_folder = '/path/ to/ controls folder'; % it could be both healthy or disease controls
folder_hs_norm = '/path/ to/ control groups for norm';

load(fullfile(folder_hs_norm, 'MRF_T1_imgset.mat'))
load(fullfile(folder_hs_norm, 'MRF_T2_imgset.mat'))
load(fullfile(folder_hs_norm, 'FAST_imgset.mat'))

VBM_flag = 0; % if need VBM, put down 1; else put down 0 (if doing subtype analysis, it is recommended to add the VBM features)

%% feature generation for patients
T = readtable(fullfile(pts_folder, 'subject_data.xlsx'), 'sheet', 'FCD_subtype', 'VariableNamingRule', 'preserve');
T_cell = table2cell(T);

disp('feature generation for patients')
subtype_flag = 0; % if not sybtype analysis, put down 0; if subtype analysis, enter 1 and calculate additional features (entropy and uniformity)
sub_all = dir(fullfile(pts_folder, 'P*'));
for order = 1:length(sub_all)
    T1_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'n_syN_T1_data_brain_Warped.nii'));
    T1_img = double(T1_nii.img);
    T2_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'n_syN_T2_data_brain_Warped.nii'));
    T2_img = double(T2_nii.img);

    ROI_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'n_FCD_ROI_bin.nii'));
    ROI_img = double(ROI_nii.img);

    CSF_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'n_syN_T1w_data_brain_Warped_seg_0.nii.gz'));
    CSF_img = double(CSF_nii.img);
    GM_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'n_syN_T1w_data_brain_Warped_seg_1.nii.gz'));
    GM_img = double(GM_nii.img);
    WM_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'n_syN_T1w_data_brain_Warped_seg_2.nii.gz'));
    WM_img = double(WM_nii.img);

    CSF_dil_img = zeros(size(CSF_img));
    SE = strel('disk', 1);
    for slice = 1:size(CSF_img, 3)
        CSF_dil_img(:, :, slice) = imdilate(CSF_img(:, :, slice), SE);
    end
    CSF_dil_img(find(T1w_img == 0)) = 1;
    
    ROI_GM_img = zeros(size(ROI_img));
    ROI_WM_img = zeros(size(ROI_img));

    ROI_GM_img(find(ROI_img + GM_img == 2)) = 1;
    ROI_WM_img(find(ROI_img + WM_img == 2)) = 1;

    ROI_img(find(CSF_dil_img)) = 0;
    ROI_GM_img(find(CSF_dil_img)) = 0;
    ROI_WM_img(find(CSF_dil_img)) = 0;
    
    [T1_GM_voxels, T1_WM_voxels, T1_norm_para] = get_MRF_value(T1_img, ROI_img, ROI_GM_img, ROI_WM_img, MRF_T1_imgset, FAST_imgset);
    [T2_GM_voxels, T2_WM_voxels, T2_norm_para] = get_MRF_value(T2_img, ROI_img, ROI_GM_img, ROI_WM_img, MRF_T2_imgset, FAST_imgset);

    [T1_fts_2d_norm, T1_fts_2d, T1_fts_3d_norm, T1_fts_3d, slice_num] = MRF_feature_gen_Indnorm(T1_img, ROI_img, ROI_GM_img, ROI_WM_img, MRF_T1_imgset, FAST_imgset, subtype_flag);
    [T2_fts_2d_norm, T2_fts_2d, T2_fts_3d_norm, T2_fts_3d, ~] = MRF_feature_gen_Indnorm(T2_img, ROI_img, ROI_GM_img, ROI_WM_img, MRF_T2_imgset, FAST_imgset, subtype_flag);

    ROI_pt_source = repmat({sub_all(order).name}, slice_num, 1);

    % MAP VBM features
    if VBM_flag == 1
        VBM_junction_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'MNI_junction_fn.nii'));
        VBM_junction_img = double(VBM_junction_nii.img);
        VBM_thickness_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'MNI_thickness_fn.nii'));
        VBM_thickness_img = double(VBM_thickness_nii.img);
        VBM_extension_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'MNI_extension_fn.nii'));
        VBM_extension_img = double(VBM_extension_nii.img);

        [VBM_fts_2d, VBM_fts_3d, ~] = VBM_feature_gen(VBM_junction_img, VBM_thickness_img, VBM_extension_img, ROI_img, ROI_GM_img, ROI_WM_img);

        norm_fts_2d = [T1_fts_2d_norm, T2_fts_2d_norm, VBM_fts_2d];
        fts_2d = [T1_fts_2d, T2_fts_2d, VBM_fts_2d];
        
        norm_fts_3d = [T1_fts_3d_norm, T2_fts_3d_norm, VBM_fts_3d];
        fts_3d = [T1_fts_3d, T2_fts_3d, VBM_fts_3d];
    else
        norm_fts_2d = [T1_fts_2d_norm, T2_fts_2d_norm];
        fts_2d = [T1_fts_2d, T2_fts_2d];
        
        norm_fts_3d = [T1_fts_3d_norm, T2_fts_3d_norm];
        fts_3d = [T1_fts_3d, T2_fts_3d];
    end

    rm_slice_num = length(find(sum(isnan(norm_fts_2d) + isinf(norm_fts_2d), 2)));

    ROI_pt_source(find(sum(isnan(norm_fts_2d) + isinf(norm_fts_2d), 2))) = [];
    fts_2d(find(sum(isnan(norm_fts_2d) + isinf(norm_fts_2d), 2)), :) = [];
    norm_fts_2d(find(sum(isnan(norm_fts_2d) + isinf(norm_fts_2d), 2)), :) = [];

    FCD_study_data(order).subject_name = sub_all(order).name;
    FCD_study_data(order).T1_GM_voxels = T1_GM_voxels;
    FCD_study_data(order).T1_WM_voxels = T1_WM_voxels;
    FCD_study_data(order).T2_GM_voxels = T2_GM_voxels;
    FCD_study_data(order).T2_WM_voxels = T2_WM_voxels;
    FCD_study_data(order).norm_fts_2d = norm_fts_2d;
    FCD_study_data(order).fts_2d = fts_2d;
    FCD_study_data(order).norm_fts_3d = norm_fts_3d;
    FCD_study_data(order).fts_3d = fts_3d;
    FCD_study_data(order).slice_num = slice_num - rm_slice_num;
        
    FCD_study_data(order).subtype = T_cell{find(strcmp(T_cell(:, 1), sub_all(order).name)), 2};
    FCD_study_data(order).outcome = T_cell{find(strcmp(T_cell(:, 1), sub_all(order).name)), 3};
    FCD_study_data(order).ROI_pt_source = ROI_pt_source;
    
    FCD_study_data(order).T1_norm_para = T1_norm_para;
    FCD_study_data(order).T2_norm_para = T2_norm_para;
    
    clear ROI_nii ROI_img ROI_GM_nii ROI_GM_img ROI_WM_nii ROI_WM_img T1_nii T1_img T2_nii T2_img 
    clear T1_fts_2d_norm T1_fts_2d T1_fts_3d_norm T1_fts_3d slice_num T2_fts_2d_norm T2_fts_2d T2_fts_3d_norm T2_fts_3d
    clear ROI_pt_source norm_fts_2d fts_2d rm_slice_num ROI_pt_source CSF_dil_img CSF_dil_nii
    clear T1_GM_voxels T1_WM_voxels T2_GM_voxels T2_WM_voxels T1_norm_para T2_norm_para
end
save(fullfile(pts_folder, 'FCD_study_patient_fts_Indnorm_th10.mat'), 'FCD_study_data')
clear FCD_study_data 

%% feature generation for controls
disp('feature generation for controls')
subtype_flag = 0;    
sub_all = dir(fullfile(controls_folder, '*'));
pts_all = dir(fullfile(pts_folder, 'P*'));
for order = 1:length(sub_all)
    T1_nii = load_untouch_nii(fullfile(controls_folder, sub_all(order).name, 'n_syN_T1_data_brain_Warped.nii'));
    T1_img = double(T1_nii.img);
    T2_nii = load_untouch_nii(fullfile(controls_folder, sub_all(order).name, 'n_syN_T2_data_brain_Warped.nii'));
    T2_img = double(T2_nii.img);
    
    CSF_nii = load_untouch_nii(fullfile(pts_folder, sub_all(order).name, 'n_syN_T1w_data_brain_Warped_seg_0.nii.gz'));
    CSF_img = double(CSF_nii.img);

    CSF_dil_img = zeros(size(CSF_img));
    SE = strel('disk', 1);
    for slice = 1:size(CSF_img, 3)
        CSF_dil_img(:, :, slice) = imdilate(CSF_img(:, :, slice), SE);
    end
    CSF_dil_img(find(T1w_img == 0)) = 1;
    
    for ROI_lbl = 1:length(pts_all)
        ROI_nii = load_untouch_nii(fullfile(pts_folder, pts_all(ROI_lbl).name, 'n_FCD_ROI_bin.nii'));
        ROI_img = double(ROI_nii.img);
        GM_nii = load_untouch_nii(fullfile(pts_folder, pts_all(ROI_lbl).name, 'n_syN_T1w_data_brain_Warped_seg_1.nii.gz'));
        GM_img = double(GM_nii.img);
        WM_nii = load_untouch_nii(fullfile(pts_folder, pts_all(ROI_lbl).name, 'n_syN_T1w_data_brain_Warped_seg_2.nii.gz'));
        WM_img = double(WM_nii.img);

        ROI_GM_img = zeros(size(ROI_img));
        ROI_WM_img = zeros(size(ROI_img));
    
        ROI_GM_img(find(ROI_img + GM_img == 2)) = 1;
        ROI_WM_img(find(ROI_img + WM_img == 2)) = 1;

        ROI_img(find(CSF_dil_img)) = 0;
        ROI_GM_img(find(CSF_dil_img)) = 0;
        ROI_WM_img(find(CSF_dil_img)) = 0;
    
        [T1_GM_voxels_pre, T1_WM_voxels_pre, T1_norm_para_pre] = get_MRF_value(T1_img, ROI_img, ROI_GM_img, ROI_WM_img, MRF_T1_imgset, FAST_imgset);
        [T2_GM_voxels_pre, T2_WM_voxels_pre, T2_norm_para_pre] = get_MRF_value(T2_img, ROI_img, ROI_GM_img, ROI_WM_img, MRF_T2_imgset, FAST_imgset);
    
        [T1_fts_2d_norm, T1_fts_2d, T1_fts_3d_norm, T1_fts_3d, slice_num] = MRF_feature_gen_Indnorm(T1_img, ROI_img, ROI_GM_img, ROI_WM_img, MRF_T1_imgset, FAST_imgset, subtype_flag);
        [T2_fts_2d_norm, T2_fts_2d, T2_fts_3d_norm, T2_fts_3d, ~] = MRF_feature_gen_Indnorm(T2_img, ROI_img, ROI_GM_img, ROI_WM_img, MRF_T2_imgset, FAST_imgset, subtype_flag);

        ROI_pt_source = repmat({pts_all(ROI_lbl).name}, slice_num, 1);
        
        if ROI_lbl == 1
            T1_fts_2d_norm_hcs = T1_fts_2d_norm;
            T2_fts_2d_norm_hcs = T2_fts_2d_norm;
            T1_fts_2d_hcs = T1_fts_2d;
            T2_fts_2d_hcs = T2_fts_2d;
            
            T1_fts_3d_norm_hcs = T1_fts_3d_norm;
            T2_fts_3d_norm_hcs = T2_fts_3d_norm;
            T1_fts_3d_hcs = T1_fts_3d;
            T2_fts_3d_hcs = T2_fts_3d;
            
            ROI_pt_source_hcs = ROI_pt_source;
            slice_num_hcs = slice_num;
            
            T1_GM_voxels = T1_GM_voxels_pre;
            T1_WM_voxels = T1_WM_voxels_pre;
            T2_GM_voxels = T2_GM_voxels_pre;
            T2_WM_voxels = T2_WM_voxels_pre;
            
            T1_vol_ROI_GM_vxls = T1_norm_para_pre.vol_ROI_GM_vxls;
            T1_vol_ROI_WM_vxls = T1_norm_para_pre.vol_ROI_WM_vxls;
            T2_vol_ROI_GM_vxls = T2_norm_para_pre.vol_ROI_GM_vxls;
            T2_vol_ROI_WM_vxls = T2_norm_para_pre.vol_ROI_WM_vxls;
        else
            T1_fts_2d_norm_hcs = [T1_fts_2d_norm_hcs; T1_fts_2d_norm];
            T2_fts_2d_norm_hcs = [T2_fts_2d_norm_hcs; T2_fts_2d_norm];
            T1_fts_2d_hcs = [T1_fts_2d_hcs; T1_fts_2d];
            T2_fts_2d_hcs = [T2_fts_2d_hcs; T2_fts_2d];
            
            T1_fts_3d_norm_hcs = [T1_fts_3d_norm_hcs; T1_fts_3d_norm];
            T2_fts_3d_norm_hcs = [T2_fts_3d_norm_hcs; T2_fts_3d_norm];
            T1_fts_3d_hcs = [T1_fts_3d_hcs; T1_fts_3d];
            T2_fts_3d_hcs = [T2_fts_3d_hcs; T2_fts_3d];
            
            ROI_pt_source_hcs = [ROI_pt_source_hcs; ROI_pt_source];
            slice_num_hcs = slice_num_hcs + slice_num;
            
            T1_GM_voxels = [T1_GM_voxels; T1_GM_voxels_pre];
            T1_WM_voxels = [T1_WM_voxels; T1_WM_voxels_pre];
            T2_GM_voxels = [T2_GM_voxels; T2_GM_voxels_pre];
            T2_WM_voxels = [T2_WM_voxels; T2_WM_voxels_pre];

            T1_vol_ROI_GM_vxls = [T1_vol_ROI_GM_vxls; T1_norm_para_pre.vol_ROI_GM_vxls];
            T1_vol_ROI_WM_vxls = [T1_vol_ROI_WM_vxls; T1_norm_para_pre.vol_ROI_WM_vxls];
            T2_vol_ROI_GM_vxls = [T2_vol_ROI_GM_vxls; T2_norm_para_pre.vol_ROI_GM_vxls];
            T2_vol_ROI_WM_vxls = [T2_vol_ROI_WM_vxls; T2_norm_para_pre.vol_ROI_WM_vxls];
        end 
        
        clear T1_GM_voxels_pre T1_WM_voxels_pre T2_GM_voxels_pre T2_WM_voxels_pre T1_norm_para_pre T2_norm_para_pre
        clear ROI_nii ROI_img ROI_GM_nii ROI_GM_img ROI_WM_nii ROI_WM_img ROI_pt_source
        clear T1_fts_2d_norm T1_fts_2d T1_fts_3d_norm T1_fts_3d T2_fts_2d_norm T2_fts_2d T2_fts_3d_norm T2_fts_3d slice_num
    end
    T1_norm_para.vol_ROI_GM_vxls = T1_vol_ROI_GM_vxls;
    T1_norm_para.vol_ROI_WM_vxls = T1_vol_ROI_WM_vxls;
    T2_norm_para.vol_ROI_GM_vxls = T2_vol_ROI_GM_vxls;
    T2_norm_para.vol_ROI_WM_vxls = T2_vol_ROI_WM_vxls;
    clear T1_vol_ROI_GM_vxls T1_vol_ROI_WM_vxls T2_vol_ROI_GM_vxls T2_vol_ROI_WM_vxls
    
    fts_2d_hcs_norm = [T1_fts_2d_norm_hcs, T2_fts_2d_norm_hcs];
    fts_2d_hcs = [T1_fts_2d_hcs, T2_fts_2d_hcs];
    
    fts_3d_hcs_norm = [T1_fts_3d_norm_hcs, T2_fts_3d_norm_hcs];
    fts_3d_hcs = [T1_fts_3d_hcs, T2_fts_3d_hcs];
    
    rm_slice_num = length(find(sum(isnan(fts_2d_hcs_norm) + isinf(fts_2d_hcs_norm), 2)));
    ROI_pt_source_hcs(find(sum(isnan(fts_2d_hcs_norm) + isinf(fts_2d_hcs_norm), 2))) = [];
    fts_2d_hcs_norm(find(sum(isnan(fts_2d_hcs_norm) + isinf(fts_2d_hcs_norm), 2)), :) = [];
    fts_2d_hcs(find(sum(isnan(fts_2d_hcs_norm) + isinf(fts_2d_hcs_norm), 2)), :) = [];
    
    FCD_study_data(order).subject_name = sub_all(order).name;
    FCD_study_data(order).T1_GM_voxels = T1_GM_voxels;
    FCD_study_data(order).T1_WM_voxels = T1_WM_voxels;
    FCD_study_data(order).T2_GM_voxels = T2_GM_voxels;
    FCD_study_data(order).T2_WM_voxels = T2_WM_voxels;
    FCD_study_data(order).norm_fts_2d = fts_2d_hcs_norm;
    FCD_study_data(order).fts_2d = fts_2d_hcs;
    FCD_study_data(order).norm_fts_3d = fts_3d_hcs_norm;
    FCD_study_data(order).fts_3d = fts_3d_hcs;
    FCD_study_data(order).slice_num = slice_num_hcs - rm_slice_num;
    FCD_study_data(order).subtype = 'Normal';
    FCD_study_data(order).ROI_pt_source = ROI_pt_source_hcs;
    FCD_study_data(order).T1_norm_para = T1_norm_para;
    FCD_study_data(order).T2_norm_para = T2_norm_para;
    
    clear T1_nii T1_img T2_nii T2_img T1_fts_hcs T2_fts_hcs slice_num_hcs fts_hcs ROI_pt_source_hcs rm_slice_num fts_hcs_norm T1_fts_norm_hcs T2_fts_norm_hcs
    clear fts_3d_hcs_norm T1_fts_3d_norm_hcs T2_fts_3d_norm_hcs fts_3d_hcs T1_fts_3d_hcs T2_fts_3d_hcs
    clear fts_2d_hcs_norm T1_fts_2d_norm_hcs T2_fts_2d_norm_hcs fts_2d_hcs T1_fts_2d_hcs T2_fts_2d_hcs
    clear rm_slice_num ROI_pt_source_hcs fts_2d_hcs_norm fts_2d_hcs CSF_dil_img CSF_dil_nii
    clear T1_GM_voxels T1_WM_voxels T2_GM_voxels T2_WM_voxels T1_norm_para T2_norm_para
end
save(fullfile(controls_folder, 'FCD_study_hcs_fts_Indnorm_th10.mat'), 'FCD_study_data')
clear FCD_study_data 
