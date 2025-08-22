clear all; close all;

toolbox_path = '/path/ to/ toolbox';

addpath(fullfile(toolbox_path, 'NIfTI_20140122'))

folder_pt = '/path/ to/ patient folder';
folder_dcs_test = '/path/ to/ disease controls folder'; 
folder_hcs_test = '/path/ to/ healthy controls folder'; 
folder_output = '/path/ to/ output folder';

%% classification (Pt vs hc)
load(fullfile(folder_hcs_test, 'FCD_study_hcs_fts_Indnorm_th10.mat'))
FCD_study_data_hcs = FCD_study_data;
clear FCD_study_data

load(fullfile(folder_pt, 'FCD_study_patient_fts_Indnorm_th10.mat'))
FCD_study_data_pts = FCD_study_data;
clear FCD_study_data

test_time = 10;

% get the list for the result evaluation
for order = 1:length(FCD_study_data_pts)
    FCD_subs_all{order, 1} = FCD_study_data_pts(order).subject_name;
    if strcmp(FCD_study_data_pts(order).subject_name(1), 'P')
        FCD_lbl_all(order, 1) = 1;
    elseif strcmp(FCD_study_data_pts(order).subject_name(1), 'V')
        FCD_lbl_all(order, 1) = 0;
    end
end

offset = length(FCD_subs_all);
for order = 1:length(FCD_study_data_hcs)
    FCD_subs_all{offset + order, 1} = FCD_study_data_hcs(order).subject_name;
    if strcmp(FCD_study_data_hcs(order).subject_name(1), 'P')
        FCD_lbl_all(offset + order, 1) = 1;
    elseif strcmp(FCD_study_data_hcs(order).subject_name(1), 'V')
        FCD_lbl_all(offset + order, 1) = 0;
    end
end

for test_order= 1:test_time
    disp(['times: ' num2str(test_order)])
%     FCD_pred_all = zeros(size(FCD_lbl_all));
    FCD_prob_all = zeros(size(FCD_lbl_all));%%
    
    % rng('default') % For reproducibility
    cvp_pts = cvpartition(length(FCD_study_data_pts), 'KFold', 5);
    cvp_hcs = cvpartition(length(FCD_study_data_hcs), 'KFold', 5);
    
    for fold_order = 1:5
        Train_pts_cvp = pts_all(training(cvp_pts, fold_order));
        Test_pts_cvp = pts_all(test(cvp_pts, fold_order));

        Train_hcs_cvp = hcs_test_all(training(cvp_hcs, fold_order));
        Test_hcs_cvp = hcs_test_all(test(cvp_hcs, fold_order));
        
        % train data generation (2D)
        for sub_order = 1:length(Train_pts_cvp)
            fts_pt_idx = find(strcmp({FCD_study_data_pts.subject_name}, Train_pts_cvp(sub_order).name));
            Train_data_summary(sub_order) = rmfield(FCD_study_data_pts(fts_pt_idx), 'outcome');
            clear fts_pt_idx
        end

        offset = length(Train_data_summary);
        for sub_order = 1:length(Train_hcs_cvp)
            fts_hc_idx = find(strcmp({FCD_study_data_hcs.subject_name}, Train_hcs_cvp(sub_order).name));
            Train_data_summary(offset + sub_order) = FCD_study_data_hcs(fts_hc_idx);
            clear fts_hc_idx
        end

        Train_data = [];
        Train_label = [];
        for sub_order = 1:length(Train_data_summary)
            Train_data = [Train_data; Train_data_summary(sub_order).norm_fts_2d];

            if strcmp(Train_data_summary(sub_order).subject_name(1), 'P')
                Train_label = [Train_label; ones(size(Train_data_summary(sub_order).norm_fts_2d, 1), 1)];
            elseif strcmp(Train_data_summary(sub_order).subject_name(1), 'V')
                Train_label = [Train_label; zeros(size(Train_data_summary(sub_order).norm_fts_2d, 1), 1)];
            end
        end
        
        % Training step; Use RUSBoost ensemble classifier to deal with the
        % imbalanced data issue
        %rng('default');
        t = templateTree('MaxNumSplits', size(Train_data, 1), 'Surrogate', 'on');
        rusTree = fitcensemble(Train_data, Train_label, 'Method', 'RUSBoost', ...
            'NumLearningCycles', 1000, 'Learners', t,'LearnRate', 0.1, 'nprint', 100);
        
        % decide the threshold from the training dataset
        [Training_pred, Training_scores] = predict(rusTree, Train_data);
        Training_prob = Training_scores(:, 2)./sum(Training_scores, 2);
        T = 0:0.01:1;
        for T_order = 1:length(T)
            Train_pred = double(Training_prob > T(T_order));
            sen(1, T_order) = length(find(Train_label + Train_pred == 2))/length(find(Train_label));
            spe(1, T_order) = length(find(Train_label + Train_pred == 0))/length(find(~Train_label));
            clear Train_pred
        end
        Youden_J = spe + sen - 1;
        threshold_2D = mean(T(find(Youden_J == max(Youden_J)))); %%%%%% important
        clear spe sen Youden_J T Training_prob Training_pred Training_scores
        threshold_2D_set(fold_order, test_order) = threshold_2D; %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % calculate the significance
        [imp, ma] = predictorImportance(rusTree);
        imp_2d(:, fold_order, test_order) = imp';
        ma_2d(:, :, fold_order, test_order) = ma;
        clear imp ma
        
        % test data generation (2D)
        for sub_order = 1:length(Test_pts_cvp)
            fts_pt_idx = find(strcmp({FCD_study_data_pts.subject_name}, Test_pts_cvp(sub_order).name));
            Test_data_summary(sub_order) = rmfield(FCD_study_data_pts(fts_pt_idx), 'outcome');
            clear fts_pt_idx
        end

        offset = length(Test_data_summary);
        for sub_order = 1:length(Test_hcs_cvp)
            fts_hc_idx = find(strcmp({FCD_study_data_hcs.subject_name}, Test_hcs_cvp(sub_order).name));
            Test_data_summary(offset + sub_order) = FCD_study_data_hcs(fts_hc_idx);
            clear fts_hc_idx
        end
        
        % evaluate 2D-level performance
        for sub_order = 1:length(Test_data_summary)
            Test_data_per = Test_data_summary(sub_order).norm_fts_2d;
            [pred, scores] = predict(rusTree, Test_data_per);
            if sub_order == 1
                pred_2d_perf = double(scores(:, 2)./sum(scores, 2) > threshold_2D); %%%%% Adaptive threshold
                prob_2d_perf = scores(:, 2)./sum(scores, 2);
                if strcmp(Test_data_summary(sub_order).subtype, 'Normal')
                    true_2d_perf = zeros(size(pred));
                else
                    true_2d_perf = ones(size(pred));
                end
            else
                pred_2d_perf = [pred_2d_perf; double(scores(:, 2)./sum(scores, 2) > threshold_2D)]; %%%%% Adaptive threshold
                prob_2d_perf = [prob_2d_perf; scores(:, 2)./sum(scores, 2)];
                if strcmp(Test_data_summary(sub_order).subtype, 'Normal')
                    true_2d_perf = [true_2d_perf; zeros(size(pred))];
                else
                    true_2d_perf = [true_2d_perf; ones(size(pred))];
                end
            end
            clear Test_data_per pred scores
        end
        
        C = confusionchart(true_2d_perf, pred_2d_perf);
        perf_2d(1:4, fold_order, test_order) = reshape(C.NormalizedValues, 4, 1);
        
        [FPR_2d_perf, TPR_2d_perf, T_2d_perf, AUC_2d_perf] = perfcurve(true_2d_perf, pred_2d_perf, 1);
        
        perf_2d(5, fold_order, test_order) = AUC_2d_perf;
        
        FPR_samp_xpin = [0:0.01:1]';
        [~, idx_uni] = unique(FPR_2d_perf);
        TPR_samp_xpin_2d(:, fold_order, test_order) = interp1(FPR_2d_perf(idx_uni), TPR_2d_perf(idx_uni), FPR_samp_xpin, 'linear');
        FPR_samp_xpin_2d(:, fold_order, test_order) = FPR_samp_xpin;
        clear idx_uni FPR_samp_xpin

        TPR_samp_ypin = [0:0.01:1]';
        [~, idx_uni] = unique(TPR_2d_perf);
        FPR_samp_ypin_2d(:, fold_order, test_order) = interp1(TPR_2d_perf(idx_uni), FPR_2d_perf(idx_uni), TPR_samp_ypin, 'linear');
        TPR_samp_ypin_2d(:, fold_order, test_order) = TPR_samp_ypin;
        clear idx_uni TPR_samp_ypin

        clear FPR_2d_perf TPR_2d_perf T_2d_perf AUC_2d_perf C

        % evaluate patient-level performance
        for sub_order = 1:length(Test_data_summary)
            Test_data_per = Test_data_summary(sub_order).norm_fts_2d;
            [~, scores] = predict(rusTree, Test_data_per);
    %         Pred_target_cluster_lbl = double(scores(:, 2) - scores(:, 1) > threshold2); 
            pred = double(scores(:, 2)./sum(scores, 2) > threshold_2D); %%%%% Adaptive threshold
            
            if strcmp(Test_data_summary(sub_order).subject_name(1), 'P')
%                 if sum(pred) >= length(pred)*threshold_2d23d
%                     FCD_pred_all(find(strcmp(FCD_subs_all, Test_data_summary(sub_order).subject_name))) = 1;
%                 else
%                     FCD_pred_all(find(strcmp(FCD_subs_all, Test_data_summary(sub_order).subject_name))) = 0;
%                 end
                FCD_prob_all(find(strcmp(FCD_subs_all, Test_data_summary(sub_order).subject_name))) = sum(pred)/length(pred);
                
            elseif strcmp(Test_data_summary(sub_order).subject_name(1), 'V')
                Test_ROI_pt_source = Test_data_summary(sub_order).ROI_pt_source;
                Test_ROI_pt_source_sets = unique(Test_ROI_pt_source);
                
                for Test_ROI_pt_source_order = 1:length(Test_ROI_pt_source_sets)
                    idx = find(strcmp(Test_ROI_pt_source, Test_ROI_pt_source_sets(Test_ROI_pt_source_order)));
                    pred_ROI_per = pred(idx);
                    
                    FCD_prob_ROI_per(Test_ROI_pt_source_order) = sum(pred_ROI_per)/length(pred_ROI_per);
                    
%                     if sum(pred_ROI_per) >= length(pred_ROI_per)*threshold_2d23d
%                         FCD_pred_ROI_per(Test_ROI_pt_source_order) = 1;
%                     else
%                         FCD_pred_ROI_per(Test_ROI_pt_source_order) = 0;
%                     end
                    clear idx pred_ROI_per
                end
                
%                 if sum(FCD_pred_ROI_per) == 0
%                     FCD_pred_all(find(strcmp(FCD_subs_all, Test_data_summary(sub_order).subject_name))) = 0;
%                 else
%                     FCD_pred_all(find(strcmp(FCD_subs_all, Test_data_summary(sub_order).subject_name))) = 1;
%                 end
                FCD_prob_all(find(strcmp(FCD_subs_all, Test_data_summary(sub_order).subject_name))) = max(FCD_prob_ROI_per);
                
                clear Test_ROI_pt_source Test_ROI_pt_source_sets FCD_pred_ROI_per FCD_prob_ROI_per
            end
            
            clear Test_data_per pred scores Test_ROI_pt_source Test_ROI_pt_source_sets
        end

        clear Train_pts_cvp Test_pts_cvp Train_hcs_cvp Test_hcs_cvp Train_data_summary 
        clear Train_data Train_label t rusTree Test_data_summary
    end
       
    %2d ROC plot   
    FPR_samp_xpin_2d_cv = mean(FPR_samp_xpin_2d(:, :, test_order), 2);
    TPR_samp_xpin_2d_mean = mean(TPR_samp_xpin_2d(:, :, test_order), 2);
    TPR_samp_xpin_2d_std = std(TPR_samp_xpin_2d(:, :, test_order), 0, 2);

    errorbar(FPR_samp_xpin_2d_cv, TPR_samp_xpin_2d_mean, TPR_samp_xpin_2d_std);
    xlim([0 1])
    ylim([0 1.1])
    title('ROC plot evaXdir');
    saveas(gcf, fullfile(folder_hcs_test, ['ROC_plot_evaXdir_PT_HC_2d_rep' num2str(test_order) '_Indnorm_th10.png']));
    close all;

    TPR_samp_ypin_2d_cv = mean(TPR_samp_ypin_2d(:, :, test_order), 2);
    FPR_samp_ypin_2d_mean = mean(FPR_samp_ypin_2d(:, :, test_order), 2);
    FPR_samp_ypin_2d_std = std(FPR_samp_ypin_2d(:, :, test_order), 0, 2);

    errorbar(FPR_samp_ypin_2d_mean, TPR_samp_ypin_2d_cv, FPR_samp_ypin_2d_std, 'horizontal');
    xlim([0 1])
    ylim([0 1.1])
    title('ROC plot evaYdir');
    saveas(gcf, fullfile(folder_hcs_test, ['ROC_plot_evaYdir_PT_HC_2d_rep' num2str(test_order) '_Indnorm_th10.png']));
    close all;
   
    clear FPR_samp_xpin_2d_cv TPR_samp_xpin_2d_mean TPR_samp_xpin_2d_std
    clear TPR_samp_ypin_2d_cv FPR_samp_ypin_2d_mean FPR_samp_ypin_2d_std
    
    % 3d ROC plot
    T = 0:0.01:1;
    for T_order = 1:length(T)
        pred = double(FCD_prob_all > T(T_order));
        sen(1, T_order) = length(find(FCD_lbl_all + pred == 2))/length(find(FCD_lbl_all));
        spe(1, T_order) = length(find(FCD_lbl_all + pred == 0))/length(find(~FCD_lbl_all));
        clear Train_pred
    end
    Youden_J = spe + sen - 1;
    threshold_3D = mean(T(find(Youden_J == max(Youden_J)))); %%%%%% important
    clear T sen spe Youden_J
    
    perf(1, test_order) = threshold_3D; % optimized threshold
    
    FCD_pred_all = double(FCD_prob_all > threshold_3D);
    perf(2, test_order) = length(find(FCD_lbl_all + FCD_pred_all == 2))/length(find(FCD_lbl_all)); % optimized sensitivity
    perf(3, test_order) = length(find(FCD_lbl_all + FCD_pred_all == 0))/length(find(~FCD_lbl_all)); % optimized specificity
    
    % confusion matrix
    C = confusionchart(FCD_lbl_all, FCD_pred_all);
    perf(4:7, test_order) = reshape(C.NormalizedValues, 4, 1); % optimized confusion matrix
    
    title(num2str(test_order));
    saveas(gcf, fullfile(folder_hcs_test, [num2str(test_order) '_Indnorm_th10.png']))
    
    [FPR, TPR, T, AUC] = perfcurve(FCD_lbl_all, FCD_prob_all, 1);
    
    FPR_samp_xpin = [0:0.01:1]';
    [~, idx_uni] = unique(FPR);
    TPR_samp_xpin(:, test_order) = interp1(FPR(idx_uni), TPR(idx_uni), FPR_samp_xpin, 'linear');
    clear idx_uni

    TPR_samp_ypin = [0:0.01:1]';
    [~, idx_uni] = unique(TPR);
    FPR_samp_ypin(:, test_order) = interp1(TPR(idx_uni), FPR(idx_uni), TPR_samp_ypin, 'linear');
    clear idx_uni
    
    clear FPR TPR T

    perf(8, test_order) = AUC;
    
    FCD_pred_all_set(:, test_order) = FCD_pred_all;
    FCD_prob_all_set(:, test_order) = FCD_prob_all;
%     FCD_AUC_all_set(:, :, test_order) = [T_sample', ROC_sen, ROC_spe];
    
    clear C FCD_pred_all cvp_pts cvp_hcs
    clear FCD_prob_all T_sample AUC ROC_sen ROC_spe
    close all
end

% summarize the 2d performance
sen_2d = squeeze(perf_2d(4, :, :)./(perf_2d(4, :, :) + perf_2d(2, :, :)));
spe_2d = squeeze(perf_2d(1, :, :)./(perf_2d(1, :, :) + perf_2d(3, :, :)));
acc_2d = squeeze((perf_2d(1, :, :) + perf_2d(4, :, :))./(perf_2d(1, :, :) + perf_2d(2, :, :) + perf_2d(3, :, :) + perf_2d(4, :, :)));
auc_2d = squeeze(perf_2d(5, :, :));
for test_order= 1:test_time
    for fold_order = 1:5 
        imp_2d_cell{fold_order, test_order} = mat2str(round(imp_2d(:, fold_order, test_order), 4));
    end
end

% a={};
% b=cellfun(@str2num, a, 'UniformOutput', false);
% for order = 1:10
%     sig_ft_per_fold = b(:, order);
%     for ft_order = 1:size(sig_ft_per_fold{1}, 1)
%         c{ft_order, order} = (sig_ft_per_fold{1}(ft_order) + sig_ft_per_fold{2}(ft_order) + sig_ft_per_fold{3}(ft_order) + sig_ft_per_fold{4}(ft_order) + sig_ft_per_fold{5}(ft_order))/5;
%     end
%     clear sig_ft_per_fold
% end

T_result = cell2table([[FCD_subs_all, mat2cell(FCD_lbl_all, ones(size(FCD_lbl_all, 1), 1), 1), ...
    mat2cell(FCD_pred_all_set, ones(size(FCD_pred_all_set, 1), 1), ones(1, size(FCD_pred_all_set, 2))), ...
    mat2cell(FCD_prob_all_set, ones(size(FCD_prob_all_set, 1), 1), ones(1, size(FCD_prob_all_set, 2)))]; ...
    [mat2cell(zeros(8, 2), ones(8, 1), ones(1, 2)), mat2cell(perf, ones(size(perf, 1), 1), ones(1, size(perf, 2))), mat2cell(zeros(8, 10), ones(8, 1), ones(1, 10))]], ...
    'VariableNames', [{'SubID', 'Label'}, strcat(repmat({'Pred'}, 1, test_time), cellfun(@num2str, mat2cell(1:test_time, 1, ones(1, test_time)), 'UniformOutput', false)), strcat(repmat({'Prob'}, 1, test_time), cellfun(@num2str, mat2cell(1:test_time, 1, ones(1, test_time)), 'UniformOutput', false))]);

% for test_order = 1:test_time
%     if test_order == 1
%         FCD_AUC_all_set_2d = FCD_AUC_all_set(:, :, test_order);
%     else
%         FCD_AUC_all_set_2d = [FCD_AUC_all_set_2d, FCD_AUC_all_set(:, :, test_order)];
%     end
% end

writetable(T_result, fullfile(folder_output, 'result_Indnorm_th10.xlsx'), 'WriteVariableNames', true, 'WriteRowNames', false, 'Sheet', 'PT_HC_result');  

FPR_samp_xpin;
TPR_samp_xpin_mean = mean(TPR_samp_xpin, 2);
TPR_samp_xpin_stderr = std(TPR_samp_xpin, 0, 2)/(10^0.5);

errorbar(FPR_samp_xpin, TPR_samp_xpin_mean, TPR_samp_xpin_stderr);
xlim([0 1])
ylim([0 1.1])
title('ROC plot evaXdir');
saveas(gcf, fullfile(folder_output, 'ROC_plot_evaXdir_PT_HC_Indnorm_th10.png'));
close all;
    
TPR_samp_ypin;
FPR_samp_ypin_mean = mean(FPR_samp_ypin, 2);
FPR_samp_ypin_stderr = std(FPR_samp_ypin, 0, 2)/(10^0.5);

errorbar(FPR_samp_ypin_mean, TPR_samp_ypin, TPR_samp_xpin_stderr, 'horizontal');
xlim([0 1])
ylim([0 1.1])
title('ROC plot evaYdir');
saveas(gcf, fullfile(folder_output, 'ROC_plot_evaYdir_PT_HC_Indnorm_th10.png'));
close all;

save(fullfile(folder_output, 'PT_HC_data_Indnorm_th10.mat'), 'TPR_samp_xpin_2d', 'FPR_samp_xpin_2d', 'FPR_samp_ypin_2d', 'TPR_samp_ypin_2d', 'perf_2d', 'imp_2d', 'ma_2d', ...
    'TPR_samp_ypin', 'FPR_samp_ypin', 'FPR_samp_xpin', 'TPR_samp_xpin', 'threshold_2D_set') 