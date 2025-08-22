function [ft_2d_norm, ft_2d, ft_3d_norm, ft_3d, slice_num, norm_para] = MRF_feature_gen_Indnorm(vol_3d, ROI_img, ROI_GM_img, ROI_WM_img, norm_img_imgset, norm_FAST_imgset, subtype_flag)
% subtype_flag = 1 means subtype analysis

% vol_per_slice = squeeze(sum(sum(ROI_img, 2), 1));
ROI_slice_idx = find(sum(sum(ROI_img, 2), 1) > 10);
% ROI_slice_idx = find(vol_per_slice > sum(ROI_img(:))/length(find(vol_per_slice)));

ROI_GM_slice_idx = find(sum(sum(ROI_GM_img, 2), 1) > 10);
ROI_WM_slice_idx = find(sum(sum(ROI_WM_img, 2), 1) > 10);

ROI_common_idx = intersect(intersect(ROI_slice_idx, ROI_GM_slice_idx), ROI_WM_slice_idx);
slice_num = length(ROI_common_idx);

if ~isempty(ROI_common_idx)

    % Assume that FAST segmentation label assign CSF as 1; GM as 2 and WM
    % as 3(usually); can be adjusted based on the actual labels in HC_norm
    % subjects
    norm_GM_lbl = 1*ones([size(norm_FAST_imgset, 4), 1]); %[3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3];
    norm_WM_lbl = 2*ones([size(norm_FAST_imgset, 4), 1]); %[1,2,1,1,2,1,2,2,1,2,1,1,2,1,2,2,2,2,2,2,2,1,2,2,2,1,2,2,1];
    norm_CSF_lbl = 3*ones([size(norm_FAST_imgset, 4), 1]); %[2,1,2,2,1,2,1,1,2,1,2,2,1,2,1,1,1,1,1,1,1,2,1,1,1,2,1,1,2];

    vol_ROI_GM_vxls = [];
    vol_ROI_WM_vxls = [];
    SE = strel("disk", 1);
    for order = 1:length(ROI_common_idx)
        ROI_GM_img_2d = ROI_GM_img(:, :, ROI_common_idx(order));
        ROI_WM_img_2d = ROI_WM_img(:, :, ROI_common_idx(order));
    %     vol_2d = vol_3d(:, :, ROI_common_idx(order));

        for volume_order = 1:size(norm_img_imgset, 4)
            norm_img = norm_img_imgset(:, :, ROI_common_idx(order), volume_order);
            norm_FAST = norm_FAST_imgset(:, :, ROI_common_idx(order), volume_order);

            norm_FAST_GM = zeros(size(norm_FAST));
            norm_FAST_WM = zeros(size(norm_FAST));
            norm_FAST_CSF = zeros(size(norm_FAST));

            norm_FAST_GM(find(norm_FAST == norm_GM_lbl(volume_order))) = 1;
            norm_FAST_WM(find(norm_FAST == norm_WM_lbl(volume_order))) = 1;
            norm_FAST_CSF(find(norm_FAST == norm_CSF_lbl(volume_order))) = 1;

            norm_FAST_CSF_dil = imdilate(norm_FAST_CSF, SE);

            ROI_GM = zeros(size(norm_FAST_GM));
            ROI_GM(find((ROI_GM_img_2d + norm_FAST_GM) == 2)) = 1;%
            ROI_GM(find(norm_FAST_CSF_dil)) = 0;

            ROI_WM = zeros(size(norm_FAST_WM));
            ROI_WM(find((ROI_WM_img_2d + norm_FAST_WM) == 2)) = 1;
            ROI_WM(find(norm_FAST_CSF_dil)) = 0;

            vol_ROI_GM_vxls = cat(1, vol_ROI_GM_vxls, norm_img(find(ROI_GM)));
            vol_ROI_WM_vxls = cat(1, vol_ROI_WM_vxls, norm_img(find(ROI_WM)));

            clear norm_img norm_FAST norm_FAST_GM norm_FAST_WM norm_FAST_CSF 
            clear norm_FAST_CSF_dil ROI_GM ROI_WM
        end
        clear ROI_GM_img_2d ROI_WM_img_2d vol_2d
    end

    vol_ROI_GM_vxls_mean = mean(vol_ROI_GM_vxls);
    vol_ROI_GM_vxls_std = std(vol_ROI_GM_vxls);

    vol_ROI_WM_vxls_mean = mean(vol_ROI_WM_vxls);
    vol_ROI_WM_vxls_std = std(vol_ROI_WM_vxls);

    norm_para.vol_ROI_GM_vxls_mean = vol_ROI_GM_vxls_mean;
    norm_para.vol_ROI_GM_vxls_std = vol_ROI_GM_vxls_std;
    norm_para.vol_ROI_WM_vxls_mean = vol_ROI_WM_vxls_mean;
    norm_para.vol_ROI_WM_vxls_std = vol_ROI_WM_vxls_std;
    
    for order = 1:length(ROI_common_idx)
        vol_2d = vol_3d(:, :, ROI_common_idx(order));

        ROI_GM_img_2d = ROI_GM_img(:, :, ROI_common_idx(order));
        ROI_WM_img_2d = ROI_WM_img(:, :, ROI_common_idx(order));

        % norm with ROI_GM
        vol_2d_norm_ROI_GM = (vol_2d - vol_ROI_GM_vxls_mean) / vol_ROI_GM_vxls_std;

        vol_2d_ROI_GM_mean(order, 1) = mean(vol_2d(find(ROI_GM_img_2d)));
        vol_2d_ROI_GM_std(order, 1) = std(vol_2d(find(ROI_GM_img_2d)));

        vol_2d_norm_ROI_GM_mean(order, 1) = mean(vol_2d_norm_ROI_GM(find(ROI_GM_img_2d)));
        vol_2d_norm_ROI_GM_std(order, 1) = std(vol_2d_norm_ROI_GM(find(ROI_GM_img_2d)));

        if subtype_flag == 1
            [vol_2d_norm_entropy, vol_2d_norm_uniformity] = EN_UN_generation(vol_2d_norm_ROI_GM(find(ROI_GM_img_2d)));
            [vol_2d_entropy, vol_2d_uniformity] = EN_UN_generation(vol_2d(find(ROI_GM_img_2d)));
            
            vol_2d_norm_ROI_GM_entropy(order, 1) = vol_2d_norm_entropy;
            vol_2d_norm_ROI_GM_uniformity(order, 1) = vol_2d_norm_uniformity;
            vol_2d_ROI_GM_entropy(order, 1) = vol_2d_entropy;
            vol_2d_ROI_GM_uniformity(order, 1) = vol_2d_uniformity;
            clear vol_2d_norm_entropy vol_2d_norm_uniformity vol_2d_entropy vol_2d_entropy
        end

        % norm with ROI_WM
        vol_2d_norm_ROI_WM = (vol_2d - vol_ROI_WM_vxls_mean) / vol_ROI_WM_vxls_std;

        vol_2d_ROI_WM_mean(order, 1) = mean(vol_2d(find(ROI_WM_img_2d)));
        vol_2d_ROI_WM_std(order, 1) = std(vol_2d(find(ROI_WM_img_2d)));

        vol_2d_norm_ROI_WM_mean(order, 1) = mean(vol_2d_norm_ROI_WM(find(ROI_WM_img_2d)));
        vol_2d_norm_ROI_WM_std(order, 1) = std(vol_2d_norm_ROI_WM(find(ROI_WM_img_2d)));

        if subtype_flag == 1
            [vol_2d_norm_entropy, vol_2d_norm_uniformity] = EN_UN_generation(vol_2d_norm_ROI_WM(find(ROI_WM_img_2d)));
            [vol_2d_entropy, vol_2d_uniformity] = EN_UN_generation(vol_2d(find(ROI_WM_img_2d)));
            
            vol_2d_norm_ROI_WM_entropy(order, 1) = vol_2d_norm_entropy;
            vol_2d_norm_ROI_WM_uniformity(order, 1) = vol_2d_norm_uniformity;
            vol_2d_ROI_WM_entropy(order, 1) = vol_2d_entropy;
            vol_2d_ROI_WM_uniformity(order, 1) = vol_2d_uniformity;
            clear vol_2d_norm_entropy vol_2d_norm_uniformity vol_2d_entropy vol_2d_entropy
        end
        clear vol_2d ROI_GM_img_2d ROI_WM_img_2d vol_2d_norm_ROI_GM vol_2d_norm_ROI_WM
    end

    slice_counter = 1;
    for order = 1:length(ROI_common_idx)
        vol_2d = vol_3d(:, :, ROI_common_idx(order));

        ROI_GM_img_2d = ROI_GM_img(:, :, ROI_common_idx(order));
        ROI_WM_img_2d = ROI_WM_img(:, :, ROI_common_idx(order));

        proc_vol_3d(:, :, slice_counter) = vol_2d;

        % norm with ROI_GM
        proc_ROI_GM_img_3d(:, :, slice_counter) = ROI_GM_img_2d;
        proc_vol_3d_norm_ROI_GM(:, :, slice_counter) = (vol_2d - vol_ROI_GM_vxls_mean) / vol_ROI_GM_vxls_std;

        % norm with ROI_WM
        proc_ROI_WM_img_3d(:, :, slice_counter) = ROI_WM_img_2d;
        proc_vol_3d_norm_ROI_WM(:, :, slice_counter) = (vol_2d - vol_ROI_WM_vxls_mean) / vol_ROI_WM_vxls_std;

        slice_counter = slice_counter + 1;

        clear vol_2d ROI_img_2d ROI_GM_img_2d ROI_WM_img_2d
    end

    proc_vol_3d_norm_ROI_GM(find(isinf(proc_vol_3d_norm_ROI_GM))) = 0;
    proc_vol_3d_norm_ROI_WM(find(isinf(proc_vol_3d_norm_ROI_WM))) = 0;

    vol_3d_ROI_GM_mean = mean(proc_vol_3d(find(proc_ROI_GM_img_3d)));
    vol_3d_ROI_GM_std = std(proc_vol_3d(find(proc_ROI_GM_img_3d)));

    vol_3d_norm_ROI_GM_mean = mean(proc_vol_3d_norm_ROI_GM(find(proc_ROI_GM_img_3d)));
    vol_3d_norm_ROI_GM_std = std(proc_vol_3d_norm_ROI_GM(find(proc_ROI_GM_img_3d)));

    if subtype_flag == 1
        [vol_3d_norm_ROI_GM_entropy, vol_3d_norm_ROI_GM_uniformity] = EN_UN_generation(proc_vol_3d_norm_ROI_GM(find(proc_ROI_GM_img_3d)));
        [vol_3d_ROI_GM_entropy, vol_3d_ROI_GM_uniformity] = EN_UN_generation(proc_vol_3d(find(proc_ROI_GM_img_3d)));
    end

    vol_3d_ROI_WM_mean = mean(proc_vol_3d(find(proc_ROI_WM_img_3d)));
    vol_3d_ROI_WM_std = std(proc_vol_3d(find(proc_ROI_WM_img_3d)));

    vol_3d_norm_ROI_WM_mean = mean(proc_vol_3d_norm_ROI_WM(find(proc_ROI_WM_img_3d)));
    vol_3d_norm_ROI_WM_std = std(proc_vol_3d_norm_ROI_WM(find(proc_ROI_WM_img_3d)));

    if subtype_flag == 1
        [vol_3d_norm_ROI_WM_entropy, vol_3d_norm_ROI_WM_uniformity] = EN_UN_generation(proc_vol_3d_norm_ROI_WM(find(proc_ROI_WM_img_3d)));
        [vol_3d_ROI_WM_entropy, vol_3d_ROI_WM_uniformity] = EN_UN_generation(proc_vol_3d(find(proc_ROI_WM_img_3d)));
    end

    if subtype_flag == 1
        ft_2d_norm = [vol_2d_norm_ROI_GM_mean, vol_2d_norm_ROI_GM_std, vol_2d_norm_ROI_GM_entropy, vol_2d_norm_ROI_GM_uniformity, ...
            ...
            vol_2d_norm_ROI_WM_mean, vol_2d_norm_ROI_WM_std, vol_2d_norm_ROI_WM_entropy, vol_2d_norm_ROI_WM_uniformity];

        ft_2d = [vol_2d_ROI_GM_mean, vol_2d_ROI_GM_std, vol_2d_ROI_GM_entropy, vol_2d_ROI_GM_uniformity, ...
            ...
            vol_2d_ROI_WM_mean, vol_2d_ROI_WM_std, vol_2d_ROI_WM_entropy, vol_2d_ROI_WM_uniformity];

        ft_3d_norm = [vol_3d_norm_ROI_GM_mean, vol_3d_norm_ROI_GM_std, vol_3d_norm_ROI_GM_entropy, vol_3d_norm_ROI_GM_uniformity, ...
            ...
            vol_3d_norm_ROI_WM_mean, vol_3d_norm_ROI_WM_std, vol_3d_norm_ROI_WM_entropy, vol_3d_norm_ROI_WM_uniformity];

        ft_3d = [vol_3d_ROI_GM_mean, vol_3d_ROI_GM_std, vol_3d_ROI_GM_entropy, vol_3d_ROI_GM_uniformity, ...
            ...
            vol_3d_ROI_WM_mean, vol_3d_ROI_WM_std, vol_3d_ROI_WM_entropy, vol_3d_ROI_WM_uniformity];
    else
        ft_2d_norm = [vol_2d_norm_ROI_GM_mean, vol_2d_norm_ROI_GM_std, ...
            ...
            vol_2d_norm_ROI_WM_mean, vol_2d_norm_ROI_WM_std];

        ft_2d = [vol_2d_ROI_GM_mean, vol_2d_ROI_GM_std, ...
            ...
            vol_2d_ROI_WM_mean, vol_2d_ROI_WM_std];

        ft_3d_norm = [vol_3d_norm_ROI_GM_mean, vol_3d_norm_ROI_GM_std, ...
            ...
            vol_3d_norm_ROI_WM_mean, vol_3d_norm_ROI_WM_std];

        ft_3d = [vol_3d_ROI_GM_mean, vol_3d_ROI_GM_std, ...
            ...
            vol_3d_ROI_WM_mean, vol_3d_ROI_WM_std];
    end
else
    ft_2d = [];
    ft_3d = [];
    ft_2d_norm = [];
    ft_3d_norm = [];
end