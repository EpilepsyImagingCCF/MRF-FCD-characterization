function [ft_2d, ft_3d, slice_num] = VBM_feature_gen(VBM_junction_img, VBM_thickness_img, VBM_extension_img, ROI_img,ROI_GM_img, ROI_WM_img)

ROI_slice_idx = find(sum(sum(ROI_img, 2), 1)>10);
% ROI_slice_idx = find(vol_per_slice > sum(ROI_img(:))/length(find(vol_per_slice)));

ROI_GM_slice_idx = find(sum(sum(ROI_GM_img, 2), 1)>10);
ROI_WM_slice_idx = find(sum(sum(ROI_WM_img, 2), 1)>10);

ROI_common_idx = intersect(intersect(ROI_slice_idx, ROI_GM_slice_idx), ROI_WM_slice_idx);
slice_num = length(ROI_common_idx);

if ~isempty(ROI_common_idx)
    
    for order = 1:length(ROI_common_idx)
        ROI_img_2d = ROI_img(:, :, ROI_common_idx(order));
        
        VBM_junction_img_2d = VBM_junction_img(:, :, ROI_common_idx(order));
        VBM_thickness_img_2d = VBM_thickness_img(:, :, ROI_common_idx(order));
        VBM_extension_img_2d = VBM_extension_img(:, :, ROI_common_idx(order));
        
        VBM_junction_img_2d_ROI_mean(order, 1) = mean(VBM_junction_img_2d(find(ROI_img_2d)));
        VBM_junction_img_2d_ROI_std(order, 1) = std(VBM_junction_img_2d(find(ROI_img_2d)));
    
        VBM_thickness_img_2d_ROI_mean(order, 1) = mean(VBM_thickness_img_2d(find(ROI_img_2d)));
        VBM_thickness_img_2d_ROI_std(order, 1) = std(VBM_thickness_img_2d(find(ROI_img_2d)));
        
        VBM_extension_img_2d_ROI_mean(order, 1) = mean(VBM_extension_img_2d(find(ROI_img_2d)));
        VBM_extension_img_2d_ROI_std(order, 1) = std(VBM_extension_img_2d(find(ROI_img_2d)));
    
        clear VBM_junction_img_2d VBM_thickness_img_2d VBM_extension_img_2d
    end
    
    slice_counter = 1;
    for order = 1:length(ROI_common_idx)
        ROI_img_3d(:, :, slice_counter) = ROI_img(:, :, ROI_common_idx(order));
        
        VBM_junction_img_3d(:, :, slice_counter) = VBM_junction_img(:, :, ROI_common_idx(order));
        VBM_thickness_img_3d(:, :, slice_counter) = VBM_thickness_img(:, :, ROI_common_idx(order));
        VBM_extension_img_3d(:, :, slice_counter) = VBM_extension_img(:, :, ROI_common_idx(order));
    
        slice_counter = slice_counter + 1;
        
    end
    
    VBM_junction_img_3d_ROI_mean = mean(VBM_junction_img_3d(find(ROI_img_3d)));
    VBM_junction_img_3d_ROI_std = std(VBM_junction_img_3d(find(ROI_img_3d)));
    
    VBM_thickness_img_3d_ROI_mean = mean(VBM_thickness_img_3d(find(ROI_img_3d)));
    VBM_thickness_img_3d_ROI_std = std(VBM_thickness_img_3d(find(ROI_img_3d)));
    
    VBM_extension_img_3d_ROI_mean = mean(VBM_extension_img_3d(find(ROI_img_3d)));
    VBM_extension_img_3d_ROI_std = std(VBM_extension_img_3d(find(ROI_img_3d)));
    
    ft_2d = [VBM_junction_img_2d_ROI_mean, VBM_junction_img_2d_ROI_std, ...
        VBM_thickness_img_2d_ROI_mean, VBM_thickness_img_2d_ROI_std, ...
        VBM_extension_img_2d_ROI_mean, VBM_extension_img_2d_ROI_std];
    
    ft_3d = [VBM_junction_img_3d_ROI_mean, VBM_junction_img_3d_ROI_std, ...
        VBM_thickness_img_3d_ROI_mean, VBM_thickness_img_3d_ROI_std, ...
        VBM_extension_img_3d_ROI_mean, VBM_extension_img_3d_ROI_std];
else
    ft_2d = [];
    ft_3d = [];

end
