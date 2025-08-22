clear all; close all;

toolbox_path = '/path/ to/ toolbox';

addpath(fullfile(toolbox_path, 'NIfTI_20140122'))

folder_pt = '/path/ to/ patient folder';
folder_dcs_test = '/path/ to/ disease controls folder'; 
folder_hcs_test = '/path/ to/ healthy controls folder'; 

%%
% For healthy controls
load(fullfile(folder_hcs_test, 'FCD_study_hcs_fts_Indnorm_th10.mat'))
FCD_study_data_hcs = FCD_study_data;
clear FCD_study_data

% For disease controls
load(fullfile(folder_dcs_test, 'FCD_study_dcs_fts_Indnorm_th10.mat'))
FCD_study_data_dcs = FCD_study_data;
clear FCD_study_data

% For subtype analysis
load(fullfile(folder_pt, 'FCD_study_patient_fts_subtype_Indnorm_th10.mat'))
FCD_study_data_pts = FCD_study_data;
clear FCD_study_data

FCD_I_count = 1;
FCD_II_count = 1;

for order = 1:length(FCD_study_data_pts)
    if strcmp(FCD_study_data_pts(order).subtype, 'FCD I')
        FCD_study_data_pts_FCD_I(FCD_I_count, 1) = FCD_study_data_pts(order);
        FCD_I_count = FCD_I_count + 1;
    elseif strcmp(FCD_study_data_pts(order).subtype, 'FCD IIA')   
        FCD_study_data_pts_FCD_II(FCD_II_count, 1) = FCD_study_data_pts(order);
        FCD_II_count = FCD_II_count + 1;
    elseif strcmp(FCD_study_data_pts(order).subtype, 'FCD IIB')   
        FCD_study_data_pts_FCD_II(FCD_II_count, 1) = FCD_study_data_pts(order);
        FCD_II_count = FCD_II_count + 1;
    end
end

% FCD_study_data_pts_MOGHE = FCD_study_data_pts_FCD_I(8:9);
% FCD_study_data_pts_mMCD = FCD_study_data_pts_FCD_I;
% FCD_study_data_pts_mMCD(8:9) = [];

FCD_IIA_count = 1;
FCD_IIB_count = 1;

for order = 1:length(FCD_study_data_pts)
    if strcmp(FCD_study_data_pts(order).subtype, 'FCD IIA')   
        FCD_study_data_pts_FCD_IIA(FCD_IIA_count, 1) = FCD_study_data_pts(order);
        FCD_IIA_count = FCD_IIA_count + 1;
    elseif strcmp(FCD_study_data_pts(order).subtype, 'FCD IIB')   
        FCD_study_data_pts_FCD_IIB(FCD_IIB_count, 1) = FCD_study_data_pts(order);
        FCD_IIB_count = FCD_IIB_count + 1;
    end
end
clear FCD_study_patient_fts

% For patients
load(fullfile(folder_pt, 'FCD_study_patient_fts_Indnorm_th10.mat'))
FCD_study_data_pts = FCD_study_data;
clear FCD_study_data

%% voxel-wise stats for each subject (individual-level) in a specific group (mean +- std)
for order = 1:length(FCD_study_data_pts)
    T1_GM_voxels_dist{order, 1} = [num2str(round(mean(FCD_study_data_pts(order).T1_GM_voxels), 2)) char(177) num2str(round(std(FCD_study_data_pts(order).T1_GM_voxels), 2))];
    T1_WM_voxels_dist{order, 1} = [num2str(round(mean(FCD_study_data_pts(order).T1_WM_voxels), 2)) char(177) num2str(round(std(FCD_study_data_pts(order).T1_WM_voxels), 2))];
    T2_GM_voxels_dist{order, 1} = [num2str(round(mean(FCD_study_data_pts(order).T2_GM_voxels), 2)) char(177) num2str(round(std(FCD_study_data_pts(order).T2_GM_voxels), 2))];
    T2_WM_voxels_dist{order, 1} = [num2str(round(mean(FCD_study_data_pts(order).T2_WM_voxels), 2)) char(177) num2str(round(std(FCD_study_data_pts(order).T2_WM_voxels), 2))];
    Subj_name{order, 1} = FCD_study_data_pts(order).subject_name;
end
T_stats_ind = [Subj_name, T1_GM_voxels_dist, T1_WM_voxels_dist, T2_GM_voxels_dist, T2_WM_voxels_dist];

%% voxel-wise stats for each group (group-level) (mean +- std)
% For patient group
for order = 1:length(FCD_study_data_pts)
    if order == 1
        T1_GM_voxels_pt = mean(rmoutliers(FCD_study_data_pts(order).T1_GM_voxels));
        T1_WM_voxels_pt = mean(rmoutliers(FCD_study_data_pts(order).T1_WM_voxels));
        T2_GM_voxels_pt = mean(rmoutliers(FCD_study_data_pts(order).T2_GM_voxels));
        T2_WM_voxels_pt = mean(rmoutliers(FCD_study_data_pts(order).T2_WM_voxels));
    else
        T1_GM_voxels_pt = cat(1, T1_GM_voxels_pt, mean(rmoutliers(FCD_study_data_pts(order).T1_GM_voxels)));
        T1_WM_voxels_pt = cat(1, T1_WM_voxels_pt, mean(rmoutliers(FCD_study_data_pts(order).T1_WM_voxels)));
        T2_GM_voxels_pt = cat(1, T2_GM_voxels_pt, mean(rmoutliers(FCD_study_data_pts(order).T2_GM_voxels)));
        T2_WM_voxels_pt = cat(1, T2_WM_voxels_pt, mean(rmoutliers(FCD_study_data_pts(order).T2_WM_voxels)));
    end
end
PT_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pt), 2)) char(177) num2str(round(std(T1_GM_voxels_pt), 2))];
PT_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pt), 2)) char(177) num2str(round(std(T1_WM_voxels_pt), 2))];
PT_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pt), 2)) char(177) num2str(round(std(T2_GM_voxels_pt), 2))];
PT_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pt), 2)) char(177) num2str(round(std(T2_WM_voxels_pt), 2))];

% For healthy controls
for order = 1:length(FCD_study_data_hcs)
    if order == 1
        T1_GM_voxels_hc = mean(rmoutliers(FCD_study_data_hcs(order).T1_GM_voxels));
        T1_WM_voxels_hc = mean(rmoutliers(FCD_study_data_hcs(order).T1_WM_voxels));
        T2_GM_voxels_hc = mean(rmoutliers(FCD_study_data_hcs(order).T2_GM_voxels));
        T2_WM_voxels_hc = mean(rmoutliers(FCD_study_data_hcs(order).T2_WM_voxels));
    else
        T1_GM_voxels_hc = cat(1, T1_GM_voxels_hc, mean(rmoutliers(FCD_study_data_hcs(order).T1_GM_voxels)));
        T1_WM_voxels_hc = cat(1, T1_WM_voxels_hc, mean(rmoutliers(FCD_study_data_hcs(order).T1_WM_voxels)));
        T2_GM_voxels_hc = cat(1, T2_GM_voxels_hc, mean(rmoutliers(FCD_study_data_hcs(order).T2_GM_voxels)));
        T2_WM_voxels_hc = cat(1, T2_WM_voxels_hc, mean(rmoutliers(FCD_study_data_hcs(order).T2_WM_voxels)));
    end
end
HC_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_hc), 2)) char(177) num2str(round(std(T1_GM_voxels_hc), 2))];
HC_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_hc), 2)) char(177) num2str(round(std(T1_WM_voxels_hc), 2))];
HC_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_hc), 2)) char(177) num2str(round(std(T2_GM_voxels_hc), 2))];
HC_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_hc), 2)) char(177) num2str(round(std(T2_WM_voxels_hc), 2))];

% For disease controls
for order = 1:length(FCD_study_data_dcs)
    if order == 1
        T1_GM_voxels_dc = mean(rmoutliers(FCD_study_data_dcs(order).T1_GM_voxels));
        T1_WM_voxels_dc = mean(rmoutliers(FCD_study_data_dcs(order).T1_WM_voxels));
        T2_GM_voxels_dc = mean(rmoutliers(FCD_study_data_dcs(order).T2_GM_voxels));
        T2_WM_voxels_dc = mean(rmoutliers(FCD_study_data_dcs(order).T2_WM_voxels));
    else
        T1_GM_voxels_dc = cat(1, T1_GM_voxels_dc, mean(rmoutliers(FCD_study_data_dcs(order).T1_GM_voxels)));
        T1_WM_voxels_dc = cat(1, T1_WM_voxels_dc, mean(rmoutliers(FCD_study_data_dcs(order).T1_WM_voxels)));
        T2_GM_voxels_dc = cat(1, T2_GM_voxels_dc, mean(rmoutliers(FCD_study_data_dcs(order).T2_GM_voxels)));
        T2_WM_voxels_dc = cat(1, T2_WM_voxels_dc, mean(rmoutliers(FCD_study_data_dcs(order).T2_WM_voxels)));
    end
end
DC_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_dc), 2)) char(177) num2str(round(std(T1_GM_voxels_dc), 2))];
DC_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_dc), 2)) char(177) num2str(round(std(T1_WM_voxels_dc), 2))];
DC_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_dc), 2)) char(177) num2str(round(std(T2_GM_voxels_dc), 2))];
DC_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_dc), 2)) char(177) num2str(round(std(T2_WM_voxels_dc), 2))];

% for order = 1:length(FCD_study_data_pts_mMCD)
%     if order == 1
%         T1_GM_voxels_pts_FCD_mMCD = mean(rmoutliers(FCD_study_data_pts_mMCD(order).T1_GM_voxels));
%         T1_WM_voxels_pts_FCD_mMCD = mean(rmoutliers(FCD_study_data_pts_mMCD(order).T1_WM_voxels));
%         T2_GM_voxels_pts_FCD_mMCD = mean(rmoutliers(FCD_study_data_pts_mMCD(order).T2_GM_voxels));
%         T2_WM_voxels_pts_FCD_mMCD = mean(rmoutliers(FCD_study_data_pts_mMCD(order).T2_WM_voxels));
%     else
%         T1_GM_voxels_pts_FCD_mMCD = cat(1, T1_GM_voxels_pts_FCD_mMCD, mean(rmoutliers(FCD_study_data_pts_mMCD(order).T1_GM_voxels)));
%         T1_WM_voxels_pts_FCD_mMCD = cat(1, T1_WM_voxels_pts_FCD_mMCD, mean(rmoutliers(FCD_study_data_pts_mMCD(order).T1_WM_voxels)));
%         T2_GM_voxels_pts_FCD_mMCD = cat(1, T2_GM_voxels_pts_FCD_mMCD, mean(rmoutliers(FCD_study_data_pts_mMCD(order).T2_GM_voxels)));
%         T2_WM_voxels_pts_FCD_mMCD = cat(1, T2_WM_voxels_pts_FCD_mMCD, mean(rmoutliers(FCD_study_data_pts_mMCD(order).T2_WM_voxels)));
%     end
% end
% pts_FCD_mMCD_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_mMCD), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_mMCD), 2))];
% pts_FCD_mMCD_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_mMCD), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_mMCD), 2))];
% pts_FCD_mMCD_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_mMCD), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_mMCD), 2))];
% pts_FCD_mMCD_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_mMCD), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_mMCD), 2))];

% for order = 1:length(FCD_study_data_pts_MOGHE)
%     if order == 1
%         T1_GM_voxels_pts_FCD_MOGHE = mean(rmoutliers(FCD_study_data_pts_MOGHE(order).T1_GM_voxels));
%         T1_WM_voxels_pts_FCD_MOGHE = mean(rmoutliers(FCD_study_data_pts_MOGHE(order).T1_WM_voxels));
%         T2_GM_voxels_pts_FCD_MOGHE = mean(rmoutliers(FCD_study_data_pts_MOGHE(order).T2_GM_voxels));
%         T2_WM_voxels_pts_FCD_MOGHE = mean(rmoutliers(FCD_study_data_pts_MOGHE(order).T2_WM_voxels));
%     else
%         T1_GM_voxels_pts_FCD_MOGHE = cat(1, T1_GM_voxels_pts_FCD_MOGHE, mean(rmoutliers(FCD_study_data_pts_MOGHE(order).T1_GM_voxels)));
%         T1_WM_voxels_pts_FCD_MOGHE = cat(1, T1_WM_voxels_pts_FCD_MOGHE, mean(rmoutliers(FCD_study_data_pts_MOGHE(order).T1_WM_voxels)));
%         T2_GM_voxels_pts_FCD_MOGHE = cat(1, T2_GM_voxels_pts_FCD_MOGHE, mean(rmoutliers(FCD_study_data_pts_MOGHE(order).T2_GM_voxels)));
%         T2_WM_voxels_pts_FCD_MOGHE = cat(1, T2_WM_voxels_pts_FCD_MOGHE, mean(rmoutliers(FCD_study_data_pts_MOGHE(order).T2_WM_voxels)));
%     end
% end
% pts_FCD_MOGHE_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_MOGHE), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_MOGHE), 2))];
% pts_FCD_MOGHE_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_MOGHE), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_MOGHE), 2))];
% pts_FCD_MOGHE_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_MOGHE), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_MOGHE), 2))];
% pts_FCD_MOGHE_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_MOGHE), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_MOGHE), 2))];

% For FCD type I group
for order = 1:length(FCD_study_data_pts_FCD_I)
    if order == 1
        T1_GM_voxels_pts_FCD_I = mean(rmoutliers(FCD_study_data_pts_FCD_I(order).T1_GM_voxels));
        T1_WM_voxels_pts_FCD_I = mean(rmoutliers(FCD_study_data_pts_FCD_I(order).T1_WM_voxels));
        T2_GM_voxels_pts_FCD_I = mean(rmoutliers(FCD_study_data_pts_FCD_I(order).T2_GM_voxels));
        T2_WM_voxels_pts_FCD_I = mean(rmoutliers(FCD_study_data_pts_FCD_I(order).T2_WM_voxels));
    else
        T1_GM_voxels_pts_FCD_I = cat(1, T1_GM_voxels_pts_FCD_I, mean(rmoutliers(FCD_study_data_pts_FCD_I(order).T1_GM_voxels)));
        T1_WM_voxels_pts_FCD_I = cat(1, T1_WM_voxels_pts_FCD_I, mean(rmoutliers(FCD_study_data_pts_FCD_I(order).T1_WM_voxels)));
        T2_GM_voxels_pts_FCD_I = cat(1, T2_GM_voxels_pts_FCD_I, mean(rmoutliers(FCD_study_data_pts_FCD_I(order).T2_GM_voxels)));
        T2_WM_voxels_pts_FCD_I = cat(1, T2_WM_voxels_pts_FCD_I, mean(rmoutliers(FCD_study_data_pts_FCD_I(order).T2_WM_voxels)));
    end
end
pts_FCD_I_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_I), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_I), 2))];
pts_FCD_I_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_I), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_I), 2))];
pts_FCD_I_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_I), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_I), 2))];
pts_FCD_I_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_I), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_I), 2))];

% For FCD type II group
for order = 1:length(FCD_study_data_pts_FCD_II)
    if order == 1
        T1_GM_voxels_pts_FCD_II = mean(rmoutliers(FCD_study_data_pts_FCD_II(order).T1_GM_voxels));
        T1_WM_voxels_pts_FCD_II = mean(rmoutliers(FCD_study_data_pts_FCD_II(order).T1_WM_voxels));
        T2_GM_voxels_pts_FCD_II = mean(rmoutliers(FCD_study_data_pts_FCD_II(order).T2_GM_voxels));
        T2_WM_voxels_pts_FCD_II = mean(rmoutliers(FCD_study_data_pts_FCD_II(order).T2_WM_voxels));
    else
        T1_GM_voxels_pts_FCD_II = cat(1, T1_GM_voxels_pts_FCD_II, mean(rmoutliers(FCD_study_data_pts_FCD_II(order).T1_GM_voxels)));
        T1_WM_voxels_pts_FCD_II = cat(1, T1_WM_voxels_pts_FCD_II, mean(rmoutliers(FCD_study_data_pts_FCD_II(order).T1_WM_voxels)));
        T2_GM_voxels_pts_FCD_II = cat(1, T2_GM_voxels_pts_FCD_II, mean(rmoutliers(FCD_study_data_pts_FCD_II(order).T2_GM_voxels)));
        T2_WM_voxels_pts_FCD_II = cat(1, T2_WM_voxels_pts_FCD_II, mean(rmoutliers(FCD_study_data_pts_FCD_II(order).T2_WM_voxels)));
    end
end
pts_FCD_II_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_II), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_II), 2))];
pts_FCD_II_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_II), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_II), 2))];
pts_FCD_II_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_II), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_II), 2))];
pts_FCD_II_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_II), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_II), 2))];

% For FCD type IIA group
for order = 1:length(FCD_study_data_pts_FCD_IIA)
    if order == 1
        T1_GM_voxels_pts_FCD_IIA = mean(rmoutliers(FCD_study_data_pts_FCD_IIA(order).T1_GM_voxels));
        T1_WM_voxels_pts_FCD_IIA = mean(rmoutliers(FCD_study_data_pts_FCD_IIA(order).T1_WM_voxels));
        T2_GM_voxels_pts_FCD_IIA = mean(rmoutliers(FCD_study_data_pts_FCD_IIA(order).T2_GM_voxels));
        T2_WM_voxels_pts_FCD_IIA = mean(rmoutliers(FCD_study_data_pts_FCD_IIA(order).T2_WM_voxels));
    else
        T1_GM_voxels_pts_FCD_IIA = cat(1, T1_GM_voxels_pts_FCD_IIA, mean(rmoutliers(FCD_study_data_pts_FCD_IIA(order).T1_GM_voxels)));
        T1_WM_voxels_pts_FCD_IIA = cat(1, T1_WM_voxels_pts_FCD_IIA, mean(rmoutliers(FCD_study_data_pts_FCD_IIA(order).T1_WM_voxels)));
        T2_GM_voxels_pts_FCD_IIA = cat(1, T2_GM_voxels_pts_FCD_IIA, mean(rmoutliers(FCD_study_data_pts_FCD_IIA(order).T2_GM_voxels)));
        T2_WM_voxels_pts_FCD_IIA = cat(1, T2_WM_voxels_pts_FCD_IIA, mean(rmoutliers(FCD_study_data_pts_FCD_IIA(order).T2_WM_voxels)));
    end
end
pts_FCD_IIA_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_IIA), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_IIA), 2))];
pts_FCD_IIA_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_IIA), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_IIA), 2))];
pts_FCD_IIA_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_IIA), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_IIA), 2))];
pts_FCD_IIA_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_IIA), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_IIA), 2))];

% For FCD type IIB group
for order = 1:length(FCD_study_data_pts_FCD_IIB)
    if order == 1
        T1_GM_voxels_pts_FCD_IIB = mean(rmoutliers(FCD_study_data_pts_FCD_IIB(order).T1_GM_voxels));
        T1_WM_voxels_pts_FCD_IIB = mean(rmoutliers(FCD_study_data_pts_FCD_IIB(order).T1_WM_voxels));
        T2_GM_voxels_pts_FCD_IIB = mean(rmoutliers(FCD_study_data_pts_FCD_IIB(order).T2_GM_voxels));
        T2_WM_voxels_pts_FCD_IIB = mean(rmoutliers(FCD_study_data_pts_FCD_IIB(order).T2_WM_voxels));
    else
        T1_GM_voxels_pts_FCD_IIB = cat(1, T1_GM_voxels_pts_FCD_IIB, mean(rmoutliers(FCD_study_data_pts_FCD_IIB(order).T1_GM_voxels)));
        T1_WM_voxels_pts_FCD_IIB = cat(1, T1_WM_voxels_pts_FCD_IIB, mean(rmoutliers(FCD_study_data_pts_FCD_IIB(order).T1_WM_voxels)));
        T2_GM_voxels_pts_FCD_IIB = cat(1, T2_GM_voxels_pts_FCD_IIB, mean(rmoutliers(FCD_study_data_pts_FCD_IIB(order).T2_GM_voxels)));
        T2_WM_voxels_pts_FCD_IIB = cat(1, T2_WM_voxels_pts_FCD_IIB, mean(rmoutliers(FCD_study_data_pts_FCD_IIB(order).T2_WM_voxels)));
    end
end
pts_FCD_IIB_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_IIB), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_IIB), 2))];
pts_FCD_IIB_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_IIB), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_IIB), 2))];
pts_FCD_IIB_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_IIB), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_IIB), 2))];
pts_FCD_IIB_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_IIB), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_IIB), 2))];

% rank sum test - PT vs HC
[T1_GM_PT_HC_p, T1_GM_PT_HC_h] = ranksum(T1_GM_voxels_pt, T1_GM_voxels_hc);
[T1_WM_PT_HC_p, T1_WM_PT_HC_h] = ranksum(T1_WM_voxels_pt, T1_WM_voxels_hc);
[T2_GM_PT_HC_p, T2_GM_PT_HC_h] = ranksum(T2_GM_voxels_pt, T2_GM_voxels_hc);
[T2_WM_PT_HC_p, T2_WM_PT_HC_h] = ranksum(T2_WM_voxels_pt, T2_WM_voxels_hc);
[T1_GM_PT_HC_p, T1_WM_PT_HC_p, T2_GM_PT_HC_p, T2_WM_PT_HC_p]

% rank sum test - PT vs DC
[T1_GM_PT_DC_p, T1_GM_PT_DC_h] = ranksum(T1_GM_voxels_pt, T1_GM_voxels_dc);
[T1_WM_PT_DC_p, T1_WM_PT_DC_h] = ranksum(T1_WM_voxels_pt, T1_WM_voxels_dc);
[T2_GM_PT_DC_p, T2_GM_PT_DC_h] = ranksum(T2_GM_voxels_pt, T2_GM_voxels_dc);
[T2_WM_PT_DC_p, T2_WM_PT_DC_h] = ranksum(T2_WM_voxels_pt, T2_WM_voxels_dc);
[T1_GM_PT_DC_p, T1_WM_PT_DC_p, T2_GM_PT_DC_p, T2_WM_PT_DC_p]

% rank sum test - FCD non-II vs II
[T1_GM_I_II_p, T1_GM_I_II_h] = ranksum(T1_GM_voxels_pts_FCD_I, T1_GM_voxels_pts_FCD_II);
[T1_WM_I_II_p, T1_WM_I_II_h] = ranksum(T1_WM_voxels_pts_FCD_I, T1_WM_voxels_pts_FCD_II);
[T2_GM_I_II_p, T2_GM_I_II_h] = ranksum(T2_GM_voxels_pts_FCD_I, T2_GM_voxels_pts_FCD_II);
[T2_WM_I_II_p, T2_WM_I_II_h] = ranksum(T2_WM_voxels_pts_FCD_I, T2_WM_voxels_pts_FCD_II);
[T1_GM_I_II_p, T1_WM_I_II_p, T2_GM_I_II_p, T2_WM_I_II_p]

% rank sum test - FCD IIA vs IIB
[T1_GM_IIA_IIB_p, T1_GM_IIA_IIB_h] = ranksum(T1_GM_voxels_pts_FCD_IIA, T1_GM_voxels_pts_FCD_IIB);
[T1_WM_IIA_IIB_p, T1_WM_IIA_IIB_h] = ranksum(T1_WM_voxels_pts_FCD_IIA, T1_WM_voxels_pts_FCD_IIB);
[T2_GM_IIA_IIB_p, T2_GM_IIA_IIB_h] = ranksum(T2_GM_voxels_pts_FCD_IIA, T2_GM_voxels_pts_FCD_IIB);
[T2_WM_IIA_IIB_p, T2_WM_IIA_IIB_h] = ranksum(T2_WM_voxels_pts_FCD_IIA, T2_WM_voxels_pts_FCD_IIB);
[T1_GM_IIA_IIB_p, T1_WM_IIA_IIB_p, T2_GM_IIA_IIB_p, T2_WM_IIA_IIB_p]

%% voxel-wise stats for each group (group-level) (mean +- std)
% Unlike previous application, this one get all the voxels from each
% subject together then calculate the mean and std

for order = 1:length(FCD_study_data_pts)
    if order == 1
        T1_GM_voxels_pt = FCD_study_data_pts(order).T1_GM_voxels;
        T1_WM_voxels_pt = FCD_study_data_pts(order).T1_WM_voxels;
        T2_GM_voxels_pt = FCD_study_data_pts(order).T2_GM_voxels;
        T2_WM_voxels_pt = FCD_study_data_pts(order).T2_WM_voxels;
    else
        T1_GM_voxels_pt = cat(1, T1_GM_voxels_pt, FCD_study_data_pts(order).T1_GM_voxels);
        T1_WM_voxels_pt = cat(1, T1_WM_voxels_pt, FCD_study_data_pts(order).T1_WM_voxels);
        T2_GM_voxels_pt = cat(1, T2_GM_voxels_pt, FCD_study_data_pts(order).T2_GM_voxels);
        T2_WM_voxels_pt = cat(1, T2_WM_voxels_pt, FCD_study_data_pts(order).T2_WM_voxels);
    end
end
PT_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pt), 2)) char(177) num2str(round(std(T1_GM_voxels_pt), 2))];
PT_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pt), 2)) char(177) num2str(round(std(T1_WM_voxels_pt), 2))];
PT_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pt), 2)) char(177) num2str(round(std(T2_GM_voxels_pt), 2))];
PT_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pt), 2)) char(177) num2str(round(std(T2_WM_voxels_pt), 2))];

for order = 1:length(FCD_study_data_hcs)
    if order == 1
        T1_GM_voxels_hc = FCD_study_data_hcs(order).T1_GM_voxels;
        T1_WM_voxels_hc = FCD_study_data_hcs(order).T1_WM_voxels;
        T2_GM_voxels_hc = FCD_study_data_hcs(order).T2_GM_voxels;
        T2_WM_voxels_hc = FCD_study_data_hcs(order).T2_WM_voxels;
    else
        T1_GM_voxels_hc = cat(1, T1_GM_voxels_hc, FCD_study_data_hcs(order).T1_GM_voxels);
        T1_WM_voxels_hc = cat(1, T1_WM_voxels_hc, FCD_study_data_hcs(order).T1_WM_voxels);
        T2_GM_voxels_hc = cat(1, T2_GM_voxels_hc, FCD_study_data_hcs(order).T2_GM_voxels);
        T2_WM_voxels_hc = cat(1, T2_WM_voxels_hc, FCD_study_data_hcs(order).T2_WM_voxels);
    end
end
HC_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_hc), 2)) char(177) num2str(round(std(T1_GM_voxels_hc), 2))];
HC_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_hc), 2)) char(177) num2str(round(std(T1_WM_voxels_hc), 2))];
HC_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_hc), 2)) char(177) num2str(round(std(T2_GM_voxels_hc), 2))];
HC_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_hc), 2)) char(177) num2str(round(std(T2_WM_voxels_hc), 2))];

for order = 1:length(FCD_study_data_dcs)
    if order == 1
        T1_GM_voxels_dc = FCD_study_data_dcs(order).T1_GM_voxels;
        T1_WM_voxels_dc = FCD_study_data_dcs(order).T1_WM_voxels;
        T2_GM_voxels_dc = FCD_study_data_dcs(order).T2_GM_voxels;
        T2_WM_voxels_dc = FCD_study_data_dcs(order).T2_WM_voxels;
    else
        T1_GM_voxels_dc = cat(1, T1_GM_voxels_dc, FCD_study_data_dcs(order).T1_GM_voxels);
        T1_WM_voxels_dc = cat(1, T1_WM_voxels_dc, FCD_study_data_dcs(order).T1_WM_voxels);
        T2_GM_voxels_dc = cat(1, T2_GM_voxels_dc, FCD_study_data_dcs(order).T2_GM_voxels);
        T2_WM_voxels_dc = cat(1, T2_WM_voxels_dc, FCD_study_data_dcs(order).T2_WM_voxels);
    end
end
DC_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_dc), 2)) char(177) num2str(round(std(T1_GM_voxels_dc), 2))];
DC_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_dc), 2)) char(177) num2str(round(std(T1_WM_voxels_dc), 2))];
DC_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_dc), 2)) char(177) num2str(round(std(T2_GM_voxels_dc), 2))];
DC_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_dc), 2)) char(177) num2str(round(std(T2_WM_voxels_dc), 2))];

for order = 1:length(FCD_study_data_pts_FCD_I)
    if order == 1
        T1_GM_voxels_pts_FCD_I = FCD_study_data_pts_FCD_I(order).T1_GM_voxels;
        T1_WM_voxels_pts_FCD_I = FCD_study_data_pts_FCD_I(order).T1_WM_voxels;
        T2_GM_voxels_pts_FCD_I = FCD_study_data_pts_FCD_I(order).T2_GM_voxels;
        T2_WM_voxels_pts_FCD_I = FCD_study_data_pts_FCD_I(order).T2_WM_voxels;
    else
        T1_GM_voxels_pts_FCD_I = cat(1, T1_GM_voxels_pts_FCD_I, FCD_study_data_pts_FCD_I(order).T1_GM_voxels);
        T1_WM_voxels_pts_FCD_I = cat(1, T1_WM_voxels_pts_FCD_I, FCD_study_data_pts_FCD_I(order).T1_WM_voxels);
        T2_GM_voxels_pts_FCD_I = cat(1, T2_GM_voxels_pts_FCD_I, FCD_study_data_pts_FCD_I(order).T2_GM_voxels);
        T2_WM_voxels_pts_FCD_I = cat(1, T2_WM_voxels_pts_FCD_I, FCD_study_data_pts_FCD_I(order).T2_WM_voxels);
    end
end
pts_FCD_I_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_I), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_I), 2))];
pts_FCD_I_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_I), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_I), 2))];
pts_FCD_I_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_I), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_I), 2))];
pts_FCD_I_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_I), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_I), 2))];

for order = 1:length(FCD_study_data_pts_FCD_II)
    if order == 1
        T1_GM_voxels_pts_FCD_II = FCD_study_data_pts_FCD_II(order).T1_GM_voxels;
        T1_WM_voxels_pts_FCD_II = FCD_study_data_pts_FCD_II(order).T1_WM_voxels;
        T2_GM_voxels_pts_FCD_II = FCD_study_data_pts_FCD_II(order).T2_GM_voxels;
        T2_WM_voxels_pts_FCD_II = FCD_study_data_pts_FCD_II(order).T2_WM_voxels;
    else
        T1_GM_voxels_pts_FCD_II = cat(1, T1_GM_voxels_pts_FCD_II, FCD_study_data_pts_FCD_II(order).T1_GM_voxels);
        T1_WM_voxels_pts_FCD_II = cat(1, T1_WM_voxels_pts_FCD_II, FCD_study_data_pts_FCD_II(order).T1_WM_voxels);
        T2_GM_voxels_pts_FCD_II = cat(1, T2_GM_voxels_pts_FCD_II, FCD_study_data_pts_FCD_II(order).T2_GM_voxels);
        T2_WM_voxels_pts_FCD_II = cat(1, T2_WM_voxels_pts_FCD_II, FCD_study_data_pts_FCD_II(order).T2_WM_voxels);
    end
end
pts_FCD_II_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_II), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_II), 2))];
pts_FCD_II_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_II), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_II), 2))];
pts_FCD_II_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_II), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_II), 2))];
pts_FCD_II_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_II), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_II), 2))];

for order = 1:length(FCD_study_data_pts_FCD_IIA)
    if order == 1
        T1_GM_voxels_pts_FCD_IIA = FCD_study_data_pts_FCD_IIA(order).T1_GM_voxels;
        T1_WM_voxels_pts_FCD_IIA = FCD_study_data_pts_FCD_IIA(order).T1_WM_voxels;
        T2_GM_voxels_pts_FCD_IIA = FCD_study_data_pts_FCD_IIA(order).T2_GM_voxels;
        T2_WM_voxels_pts_FCD_IIA = FCD_study_data_pts_FCD_IIA(order).T2_WM_voxels;
    else
        T1_GM_voxels_pts_FCD_IIA = cat(1, T1_GM_voxels_pts_FCD_IIA, FCD_study_data_pts_FCD_IIA(order).T1_GM_voxels);
        T1_WM_voxels_pts_FCD_IIA = cat(1, T1_WM_voxels_pts_FCD_IIA, FCD_study_data_pts_FCD_IIA(order).T1_WM_voxels);
        T2_GM_voxels_pts_FCD_IIA = cat(1, T2_GM_voxels_pts_FCD_IIA, FCD_study_data_pts_FCD_IIA(order).T2_GM_voxels);
        T2_WM_voxels_pts_FCD_IIA = cat(1, T2_WM_voxels_pts_FCD_IIA, FCD_study_data_pts_FCD_IIA(order).T2_WM_voxels);
    end
end
pts_FCD_IIA_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_IIA), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_IIA), 2))];
pts_FCD_IIA_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_IIA), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_IIA), 2))];
pts_FCD_IIA_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_IIA), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_IIA), 2))];
pts_FCD_IIA_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_IIA), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_IIA), 2))];

for order = 1:length(FCD_study_data_pts_FCD_IIB)
    if order == 1
        T1_GM_voxels_pts_FCD_IIB = FCD_study_data_pts_FCD_IIB(order).T1_GM_voxels;
        T1_WM_voxels_pts_FCD_IIB = FCD_study_data_pts_FCD_IIB(order).T1_WM_voxels;
        T2_GM_voxels_pts_FCD_IIB = FCD_study_data_pts_FCD_IIB(order).T2_GM_voxels;
        T2_WM_voxels_pts_FCD_IIB = FCD_study_data_pts_FCD_IIB(order).T2_WM_voxels;
    else
        T1_GM_voxels_pts_FCD_IIB = cat(1, T1_GM_voxels_pts_FCD_IIB, FCD_study_data_pts_FCD_IIB(order).T1_GM_voxels);
        T1_WM_voxels_pts_FCD_IIB = cat(1, T1_WM_voxels_pts_FCD_IIB, FCD_study_data_pts_FCD_IIB(order).T1_WM_voxels);
        T2_GM_voxels_pts_FCD_IIB = cat(1, T2_GM_voxels_pts_FCD_IIB, FCD_study_data_pts_FCD_IIB(order).T2_GM_voxels);
        T2_WM_voxels_pts_FCD_IIB = cat(1, T2_WM_voxels_pts_FCD_IIB, FCD_study_data_pts_FCD_IIB(order).T2_WM_voxels);
    end
end
pts_FCD_IIB_voxels_grp_dist{1, 1} = [num2str(round(mean(T1_GM_voxels_pts_FCD_IIB), 2)) char(177) num2str(round(std(T1_GM_voxels_pts_FCD_IIB), 2))];
pts_FCD_IIB_voxels_grp_dist{1, 2} = [num2str(round(mean(T1_WM_voxels_pts_FCD_IIB), 2)) char(177) num2str(round(std(T1_WM_voxels_pts_FCD_IIB), 2))];
pts_FCD_IIB_voxels_grp_dist{1, 3} = [num2str(round(mean(T2_GM_voxels_pts_FCD_IIB), 2)) char(177) num2str(round(std(T2_GM_voxels_pts_FCD_IIB), 2))];
pts_FCD_IIB_voxels_grp_dist{1, 4} = [num2str(round(mean(T2_WM_voxels_pts_FCD_IIB), 2)) char(177) num2str(round(std(T2_WM_voxels_pts_FCD_IIB), 2))];

% student ttest - PT vs HC
[T1_GM_PT_HC_h, T1_GM_PT_HC_p] = ttest2(T1_GM_voxels_pt, T1_GM_voxels_hc);
[T1_WM_PT_HC_h, T1_WM_PT_HC_p] = ttest2(T1_WM_voxels_pt, T1_WM_voxels_hc);
[T2_GM_PT_HC_h, T2_GM_PT_HC_p] = ttest2(T2_GM_voxels_pt, T2_GM_voxels_hc);
[T2_WM_PT_HC_h, T2_WM_PT_HC_p] = ttest2(T2_WM_voxels_pt, T2_WM_voxels_hc);
[T1_GM_PT_HC_p, T1_WM_PT_HC_p, T2_GM_PT_HC_p, T2_WM_PT_HC_p]

% student ttest - PT vs DC
[T1_GM_PT_DC_h, T1_GM_PT_DC_p] = ttest2(T1_GM_voxels_pt, T1_GM_voxels_dc);
[T1_WM_PT_DC_h, T1_WM_PT_DC_p] = ttest2(T1_WM_voxels_pt, T1_WM_voxels_dc);
[T2_GM_PT_DC_h, T2_GM_PT_DC_p] = ttest2(T2_GM_voxels_pt, T2_GM_voxels_dc);
[T2_WM_PT_DC_h, T2_WM_PT_DC_p] = ttest2(T2_WM_voxels_pt, T2_WM_voxels_dc);
[T1_GM_PT_DC_p, T1_WM_PT_DC_p, T2_GM_PT_DC_p, T2_WM_PT_DC_p]

% student ttest - FCD I vs II
[T1_GM_I_II_h, T1_GM_I_II_p] = ttest2(T1_GM_voxels_pts_FCD_I, T1_GM_voxels_pts_FCD_II);
[T1_WM_I_II_h, T1_WM_I_II_p] = ttest2(T1_WM_voxels_pts_FCD_I, T1_WM_voxels_pts_FCD_II);
[T2_GM_I_II_h, T2_GM_I_II_p] = ttest2(T2_GM_voxels_pts_FCD_I, T2_GM_voxels_pts_FCD_II);
[T2_WM_I_II_h, T2_WM_I_II_p] = ttest2(T2_WM_voxels_pts_FCD_I, T2_WM_voxels_pts_FCD_II);
[T1_GM_I_II_p, T1_WM_I_II_p, T2_GM_I_II_p, T2_WM_I_II_p]

% student ttest - FCD IIA vs IIB
[T1_GM_IIA_IIB_h, T1_GM_IIA_IIB_p] = ttest2(T1_GM_voxels_pts_FCD_IIA, T1_GM_voxels_pts_FCD_IIB);
[T1_WM_IIA_IIB_h, T1_WM_IIA_IIB_p] = ttest2(T1_WM_voxels_pts_FCD_IIA, T1_WM_voxels_pts_FCD_IIB);
[T2_GM_IIA_IIB_h, T2_GM_IIA_IIB_p] = ttest2(T2_GM_voxels_pts_FCD_IIA, T2_GM_voxels_pts_FCD_IIB);
[T2_WM_IIA_IIB_h, T2_WM_IIA_IIB_p] = ttest2(T2_WM_voxels_pts_FCD_IIA, T2_WM_voxels_pts_FCD_IIB);
[T1_GM_IIA_IIB_p, T1_WM_IIA_IIB_p, T2_GM_IIA_IIB_p, T2_WM_IIA_IIB_p]

%% FCD patients vs hc vs dc
for order = 1:length(FCD_study_data_pts)
    if order == 1
        pts_norm_fts_set = FCD_study_data_pts(order).norm_fts_2d;
        pts_fts_set = FCD_study_data_pts(order).fts_2d;
    else
        pts_norm_fts_set = cat(1, pts_norm_fts_set, FCD_study_data_pts(order).norm_fts_2d);
        pts_fts_set = cat(1, pts_fts_set, FCD_study_data_pts(order).fts_2d);
    end
end
pts_fts_set_mean = mean(pts_fts_set, 1);
pts_fts_set_std = std(pts_fts_set);

pts_norm_fts_set_mean = mean(pts_norm_fts_set, 1);
pts_norm_fts_set_std = std(pts_norm_fts_set);
    
for order = 1:length(pts_norm_fts_set_mean)
    pts_fts_set_text{order, 1} = [num2str(round(pts_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(pts_fts_set_std(order), 2)) '|' ...
        '(' num2str(round(pts_norm_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(pts_norm_fts_set_std(order), 2)) ')'];
end

%
for order = 1:length(FCD_study_data_hcs)
    if order == 1
        hcs_norm_fts_set = FCD_study_data_hcs(order).norm_fts_2d;
        hcs_fts_set = FCD_study_data_hcs(order).fts_2d;
    else
        hcs_norm_fts_set = cat(1, hcs_norm_fts_set, FCD_study_data_hcs(order).norm_fts_2d);
        hcs_fts_set = cat(1, hcs_fts_set, FCD_study_data_hcs(order).fts_2d);
    end
end

hcs_fts_set_mean = mean(hcs_fts_set, 1);
hcs_fts_set_std = std(hcs_fts_set);

hcs_norm_fts_set_mean = mean(hcs_norm_fts_set, 1);
hcs_norm_fts_set_std = std(hcs_norm_fts_set);
    
for order = 1:length(hcs_norm_fts_set_mean)
    hcs_fts_set_text{order, 1} = [num2str(round(hcs_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(hcs_fts_set_std(order), 2)) '|' ...
        '(' num2str(round(hcs_norm_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(hcs_norm_fts_set_std(order), 2)) ')'];
end

%
for order = 1:length(FCD_study_data_dcs)
    if order == 1
        dcs_norm_fts_set = FCD_study_data_dcs(order).norm_fts_2d;
        dcs_fts_set = FCD_study_data_dcs(order).fts_2d;
    else
        dcs_norm_fts_set = cat(1, dcs_norm_fts_set, FCD_study_data_dcs(order).norm_fts_2d);
        dcs_fts_set = cat(1, dcs_fts_set, FCD_study_data_dcs(order).fts_2d);
    end
end

dcs_fts_set_mean = mean(dcs_fts_set, 1);
dcs_fts_set_std = std(dcs_fts_set);

dcs_norm_fts_set_mean = mean(dcs_norm_fts_set, 1);
dcs_norm_fts_set_std = std(dcs_norm_fts_set);
    
for order = 1:length(dcs_norm_fts_set_mean)
    dcs_fts_set_text{order, 1} = [num2str(round(dcs_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(dcs_fts_set_std(order), 2)) '|' ...
        '(' num2str(round(dcs_norm_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(dcs_norm_fts_set_std(order), 2)) ')'];
end

%PT vs HC
for order = 1:length(pts_norm_fts_set_mean)
    [h, p] = ttest2(pts_fts_set(:, order), hcs_fts_set(:, order));
    [h_norm, p_norm] = ttest2(pts_norm_fts_set(:, order), hcs_norm_fts_set(:, order));
    
    p_set_PT_HC_text{order, 1} = [num2str(p) '|'  '(' num2str(p_norm) ')'];
    
%     p_set_PT_HC{order, 1} = p;
%     p_set_PT_HC{order, 2} = p_norm;
    
    clear p p_norm h h_norm
end

%PT vs DC
for order = 1:length(pts_norm_fts_set_mean)
    [h, p] = ttest2(pts_fts_set(:, order), dcs_fts_set(:, order));
    [h_norm, p_norm] = ttest2(pts_norm_fts_set(:, order), dcs_norm_fts_set(:, order));
    
    p_set_PT_DC_text{order, 1} = [num2str(p) '|'  '(' num2str(p_norm) ')'];
    
%     p_set_PT_HC{order, 1} = p;
%     p_set_PT_HC{order, 2} = p_norm;
    
    clear p p_norm h h_norm
end

%% FCD IIA vs IIB vs I
% IIA
for order = 1:length(FCD_study_data_pts_FCD_IIA)
    if order == 1
        IIA_norm_fts_set = FCD_study_data_pts_FCD_IIA(order).norm_fts_2d;
        IIA_fts_set = FCD_study_data_pts_FCD_IIA(order).fts_2d;
    else
        IIA_norm_fts_set = cat(1, IIA_norm_fts_set, FCD_study_data_pts_FCD_IIA(order).norm_fts_2d);
        IIA_fts_set = cat(1,IIA_fts_set, FCD_study_data_pts_FCD_IIA(order).fts_2d);
    end
end
IIA_fts_set(:, 17:22) = IIA_norm_fts_set(:, 17:22);

IIA_fts_set_mean = mean(IIA_fts_set, 1);
IIA_fts_set_std = std(IIA_fts_set);

IIA_norm_fts_set_mean = mean(IIA_norm_fts_set, 1);
IIA_norm_fts_set_std = std(IIA_norm_fts_set);
    
for order = 1:length(IIA_norm_fts_set_mean)
    IIA_fts_set_text{order, 1} = [num2str(round(IIA_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(IIA_fts_set_std(order), 2)) '|' ...
        '(' num2str(round(IIA_norm_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(IIA_norm_fts_set_std(order), 2)) ')'];
end

% IIB
for order = 1:length(FCD_study_data_pts_FCD_IIB)
    if order == 1
        IIB_norm_fts_set = FCD_study_data_pts_FCD_IIB(order).norm_fts_3d;
        IIB_fts_set = FCD_study_data_pts_FCD_IIB(order).fts_3d;
    else
        IIB_norm_fts_set = cat(1, IIB_norm_fts_set, FCD_study_data_pts_FCD_IIB(order).norm_fts_3d);
        IIB_fts_set = cat(1, IIB_fts_set, FCD_study_data_pts_FCD_IIB(order).fts_3d);
    end
end
IIB_fts_set(:, 17:22) = IIB_norm_fts_set(:, 17:22);

IIB_fts_set_mean = mean(IIB_fts_set, 1);
IIB_fts_set_std = std(IIB_fts_set);

IIB_norm_fts_set_mean = mean(IIB_norm_fts_set, 1);
IIB_norm_fts_set_std = std(IIB_norm_fts_set);
    
for order = 1:length(IIB_norm_fts_set_mean)
    IIB_fts_set_text{order, 1} = [num2str(round(IIB_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(IIB_fts_set_std(order), 2)) '|' ...
        '(' num2str(round(IIB_norm_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(IIB_norm_fts_set_std(order), 2)) ')'];
end

% I
for order = 1:length(FCD_study_data_pts_FCD_I)
    if order == 1
        I_norm_fts_set = FCD_study_data_pts_FCD_I(order).norm_fts_3d;
        I_fts_set = FCD_study_data_pts_FCD_I(order).fts_3d;
    else
        I_norm_fts_set = cat(1, I_norm_fts_set, FCD_study_data_pts_FCD_I(order).norm_fts_3d);
        I_fts_set = cat(1, I_fts_set, FCD_study_data_pts_FCD_I(order).fts_3d);
    end
end
I_fts_set(:, 17:22) = I_norm_fts_set(:, 17:22);

I_fts_set_mean = mean(I_fts_set, 1);
I_fts_set_std = std(I_fts_set);

I_norm_fts_set_mean = mean(I_norm_fts_set, 1);
I_norm_fts_set_std = std(I_norm_fts_set);
    
for order = 1:length(I_norm_fts_set_mean)
    I_fts_set_text{order, 1} = [num2str(round(I_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(I_fts_set_std(order), 2)) '|' ...
        '(' num2str(round(I_norm_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(I_norm_fts_set_std(order), 2)) ')'];
end

% II
for order = 1:length(FCD_study_data_pts_FCD_II)
    if order == 1
        II_norm_fts_set = FCD_study_data_pts_FCD_II(order).norm_fts_3d;
        II_fts_set = FCD_study_data_pts_FCD_II(order).fts_3d;
    else
        II_norm_fts_set = cat(1, II_norm_fts_set, FCD_study_data_pts_FCD_II(order).norm_fts_3d);
        II_fts_set = cat(1, II_fts_set, FCD_study_data_pts_FCD_II(order).fts_3d);
    end
end
II_fts_set(:, 17:22) = II_norm_fts_set(:, 17:22);

II_fts_set_mean = mean(II_fts_set, 1);
II_fts_set_std = std(II_fts_set);

II_norm_fts_set_mean = mean(II_norm_fts_set, 1);
II_norm_fts_set_std = std(II_norm_fts_set);
    
for order = 1:length(II_norm_fts_set_mean)
    II_fts_set_text{order, 1} = [num2str(round(II_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(II_fts_set_std(order), 2)) '|' ...
        '(' num2str(round(II_norm_fts_set_mean(order), 2)) ' ' char(177) ' ' num2str(round(II_norm_fts_set_std(order), 2)) ')'];
end

% II vs I
for order = 1:length(II_norm_fts_set_mean)
    [h, p] = ttest2(II_fts_set(:, order), I_fts_set(:, order));
    [h_norm, p_norm] = ttest2(II_norm_fts_set(:, order), I_norm_fts_set(:, order));
    
    p_set_II_I_text{order, 1} = [num2str(p) '|'  '(' num2str(p_norm) ')'];
    
%     p_set_PT_HC{order, 1} = p;
%     p_set_PT_HC{order, 2} = p_norm;
    
    clear p p_norm h h_norm
end

% IIb vs IIa
for order = 1:length(IIA_norm_fts_set_mean)
    [h, p] = ttest2(IIA_fts_set(:, order), IIB_fts_set(:, order));
    [h_norm, p_norm] = ttest2(IIA_norm_fts_set(:, order), IIB_norm_fts_set(:, order));
    
    p_set_IIA_IIB_text{order, 1} = [num2str(p) '|'  '(' num2str(p_norm) ')'];
    
%     p_set_PT_HC{order, 1} = p;
%     p_set_PT_HC{order, 2} = p_norm;
    
    clear p p_norm h h_norm
end