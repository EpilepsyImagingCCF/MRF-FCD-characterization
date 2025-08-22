clear all; close all;

toolbox_path = '/path/ to/ toolbox';

addpath(fullfile(toolbox_path, 'NIfTI_20140122'))

HCs_folder = '/path/ to/ control subjects';

HCs_list = dir(fullfile(HCs_folder, 'V*'));
threshold = 500;

for order = 1:length(HCs_list)
    disp([num2str(order) '/'  num2str(length(HCs_list))]);
    T1_filename = dir(fullfile(HCs_list(order).folder, HCs_list(order).name, 'n_syN_T1_*_Warped.nii'));
    T2_filename = dir(fullfile(HCs_list(order).folder, HCs_list(order).name, 'n_syN_T2_*_Warped.nii'));
    FAST_filename = dir(fullfile(HCs_list(order).folder, HCs_list(order).name, 'n_syN_*_pveseg.nii.gz'));

    MRF_T1_nii = load_untouch_nii(fullfile(T1_filename.folder, T1_filename.name));
    MRF_T2_nii = load_untouch_nii(fullfile(T2_filename.folder, T2_filename.name));
    FAST_nii = load_untouch_nii(fullfile(FAST_filename.folder, FAST_filename.name));

    MRF_T1_img = MRF_T1_nii.img;
    MRF_T2_img = MRF_T2_nii.img;
    FAST_img = FAST_nii.img;

    % intensity truncation
    MRF_T2_img(find(MRF_T2_img > threshold)) = threshold;

    if order == 1
        MRF_T1_imgset = MRF_T1_img;
        MRF_T2_imgset = MRF_T2_img;
        FAST_imgset = FAST_img;
    else
        MRF_T1_imgset = cat(4, MRF_T1_imgset, MRF_T1_img);
        MRF_T2_imgset = cat(4, MRF_T2_imgset, MRF_T2_img);
        FAST_imgset = cat(4, FAST_imgset, FAST_img);
    end

    clear T1_filename T2_filename MRF_T1_img MRF_T2_img FAST_img FAST_nii FAST_filename
end

save(fullfile(HCs_folder, 'MRF_T1_imgset.mat'), 'MRF_T1_imgset', '-v7.3')
save(fullfile(HCs_folder, 'MRF_T2_imgset.mat'), 'MRF_T2_imgset', '-v7.3')
save(fullfile(HCs_folder, 'FAST_imgset.mat'), 'FAST_imgset', '-v7.3')

%%
% HCs_MRF_T1_pveseg = load_untouch_nii(fullfile(HCs_folder, 'MRF_T1_imgset_mean_pveseg.nii.gz'));
% HCs_MRF_T1_pveseg_img = HCs_MRF_T1_pveseg.img;
% 
% HCs_MRF_T1_pveseg_CSF = zeros(size(HCs_MRF_T1_pveseg_img));
% HCs_MRF_T1_pveseg_GM = zeros(size(HCs_MRF_T1_pveseg_img));
% HCs_MRF_T1_pveseg_WM = zeros(size(HCs_MRF_T1_pveseg_img));
% 
% HCs_MRF_T1_pveseg_CSF(find(HCs_MRF_T1_pveseg_img == 1)) = 1;
% HCs_MRF_T1_pveseg_GM(find(HCs_MRF_T1_pveseg_img == 3)) = 1;
% HCs_MRF_T1_pveseg_WM(find(HCs_MRF_T1_pveseg_img == 2)) = 1;
% 
% SE = strel("sphere", 1);
% HCs_MRF_T1_pveseg_CSF_dil = imdilate(HCs_MRF_T1_pveseg_CSF, SE);
% HCs_MRF_T1_pveseg_GM(find(HCs_MRF_T1_pveseg_CSF_dil)) = 0;
% HCs_MRF_T1_pveseg_WM(find(HCs_MRF_T1_pveseg_CSF_dil)) = 0;
% 
% HCs_MRF_T1_pveseg.img = HCs_MRF_T1_pveseg_GM;
% save_untouch_nii(HCs_MRF_T1_pveseg, fullfile(HCs_folder, 'MRF_T1_imgset_mean_GM_dil.nii'))
% HCs_MRF_T1_pveseg.img = HCs_MRF_T1_pveseg_WM;
% save_untouch_nii(HCs_MRF_T1_pveseg, fullfile(HCs_folder, 'MRF_T1_imgset_mean_WM_dil.nii'))
% 
% % T2
% HCs_MRF_T2_pveseg = load_untouch_nii(fullfile(HCs_folder, 'MRF_T2_imgset_mean_pveseg.nii.gz'));
% HCs_MRF_T2_pveseg_img = HCs_MRF_T2_pveseg.img;
% 
% HCs_MRF_T2_pveseg_CSF = zeros(size(HCs_MRF_T2_pveseg_img));
% HCs_MRF_T2_pveseg_GM = zeros(size(HCs_MRF_T2_pveseg_img));
% HCs_MRF_T2_pveseg_WM = zeros(size(HCs_MRF_T2_pveseg_img));
% 
% HCs_MRF_T2_pveseg_CSF(find(HCs_MRF_T1_pveseg_img == 1)) = 1;
% HCs_MRF_T2_pveseg_GM(find(HCs_MRF_T1_pveseg_img == 3)) = 1;
% HCs_MRF_T2_pveseg_WM(find(HCs_MRF_T1_pveseg_img == 2)) = 1;
% 
% HCs_MRF_T2_pveseg_CSF_dil = imdilate(HCs_MRF_T2_pveseg_CSF, SE);
% HCs_MRF_T2_pveseg_GM(find(HCs_MRF_T2_pveseg_CSF_dil)) = 0;
% HCs_MRF_T2_pveseg_WM(find(HCs_MRF_T2_pveseg_CSF_dil)) = 0;
% 
% HCs_MRF_T2_pveseg.img = HCs_MRF_T2_pveseg_GM;
% save_untouch_nii(HCs_MRF_T2_pveseg, fullfile(HCs_folder, 'MRF_T2_imgset_mean_GM_dil.nii'))
% HCs_MRF_T2_pveseg.img = HCs_MRF_T2_pveseg_WM;
% save_untouch_nii(HCs_MRF_T2_pveseg, fullfile(HCs_folder, 'MRF_T2_imgset_mean_WM_dil.nii'))