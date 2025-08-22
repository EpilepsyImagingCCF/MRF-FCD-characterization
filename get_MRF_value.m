function [GM_voxels, WM_voxels, norm_para] = get_MRF_value(vol_3d, ROI_img, ROI_GM_img, ROI_WM_img, norm_img_imgset, norm_FAST_imgset)
% subtype_flag = 1 means subtype analysis
ROI_slice_idx = find(sum(sum(ROI_img, 2), 1) > 10);

ROI_GM_slice_idx = find(sum(sum(ROI_GM_img, 2), 1) > 10);
ROI_WM_slice_idx = find(sum(sum(ROI_WM_img, 2), 1) > 10);

ROI_common_idx = intersect(intersect(ROI_slice_idx, ROI_GM_slice_idx), ROI_WM_slice_idx);
slice_num = length(ROI_common_idx);

GM_voxels = [];
WM_voxels = [];

if ~isempty(ROI_common_idx)
    norm_img_imgset(:, :, :, 23) = [];
    norm_FAST_imgset(:, :, :, 23) = [];

    % Assume that FAST segmentation label assign CSF as 1; GM as 2 and WM
    % as 3(usually); can be adjusted based on the actual labels in HC_norm
    % subjects
    norm_CSF_lbl = 1*ones([size(norm_FAST_imgset, 4), 1]); %[2,1,2,2,1,2,1,1,2,1,2,2,1,2,1,1,1,1,1,1,1,2,1,1,1,2,1,1,2];
    norm_GM_lbl = 2*ones([size(norm_FAST_imgset, 4), 1]); %[3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3];
    norm_WM_lbl = 3*ones([size(norm_FAST_imgset, 4), 1]); %[1,2,1,1,2,1,2,2,1,2,1,1,2,1,2,2,2,2,2,2,2,1,2,2,2,1,2,2,1];
    
    vol_ROI_GM_vxls = [];
    vol_ROI_WM_vxls = [];
    SE = strel("disk", 1);
    for order = 1:length(ROI_common_idx)
        vol_2d = vol_3d(:, :, ROI_common_idx(order));

        ROI_GM_img_2d = ROI_GM_img(:, :, ROI_common_idx(order));
        ROI_WM_img_2d = ROI_WM_img(:, :, ROI_common_idx(order));

        GM_voxels = [GM_voxels; vol_2d(find(ROI_GM_img_2d))];
        WM_voxels = [WM_voxels; vol_2d(find(ROI_WM_img_2d))];

        % norm para
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
        clear vol_2d ROI_GM_img_2d ROI_WM_img_2d vol_2d_norm_ROI_GM vol_2d_norm_ROI_WM
    end

    norm_para.vol_ROI_GM_vxls = vol_ROI_GM_vxls;
    norm_para.vol_ROI_WM_vxls = vol_ROI_WM_vxls;

else
    GM_voxels = [];
    WM_voxels = [];
end