function dots = detect_dots(imdata, seg_im, num_channels, thresholds)
     MEDIAN_FILTER_SIZE = 5;
     
      % TO DO: make mask independent of binning
     if size(imdata,1) > 512 && size(imdata,2) > 512
         MEDIAN_FILTER_SIZE = 10;
     end
     
    % Only work on segmented cells to speed up processing
    [seg_im,y_ind,x_ind]= autocrop(seg_im);
    dots = struct; DOT_AREA = 4; centroids = cell(num_channels, 1);
    
    for k = 1:num_channels
        im_stack = imdata{k}(y_ind,x_ind,:);
        for p = 1:size(im_stack, 3)
            bg_sub_imstack(:, :, p) = im_stack(:, :, p) - medfilt2(im_stack(:, :, p), [MEDIAN_FILTER_SIZE MEDIAN_FILTER_SIZE]); %logMask(im_stack(:, :, p)).*uint16(logical(seg_im));
        end
        log_filt_imstack = bg_sub_imstack; %imgaussfilt3(bg_sub_imstack, 0.8); %0*im_stack;
        pot_dots_im = imregionalmax(log_filt_imstack);
                
        thresh_stack = bsxfun(@times, bg_sub_imstack, uint16(seg_im)); 
        thresh_stack(thresh_stack < thresholds(k)) = 0;
        
        % Hot pixel removal
        if size(im_stack, 3)>8
            MAX_Z_ABOVE_THRESH = round(size(im_stack, 3)/3);
            potential_hot_pixel = sum(logical(thresh_stack),3) > MAX_Z_ABOVE_THRESH;
            
            if sum(potential_hot_pixel(:))
                thresh_stack(repmat(potential_hot_pixel,1,1,size(im_stack, 3)))=0;
            end
        end
        labeled_stack = bwlabeln(thresh_stack, 26); 
        stats = regionprops(labeled_stack);
        all_dots = 1:length(stats);
        maxima_dots = unique(nonzeros(labeled_stack.*pot_dots_im));
        big_dots = all_dots([stats.Area] >= DOT_AREA);
        idx_to_keep = intersect(maxima_dots, big_dots);
        
%         conn_list = bwconncomp(thresh_stack, 26);        
%         stats = regionprops(conn_list);
%         idx_to_keep = [stats.Area] >= DOT_AREA;
                
        good_dots_stats = stats(idx_to_keep);
        dots(k).properties = good_dots_stats;
        dots(k).counts = length(good_dots_stats);  
        if dots(k).counts
                
            if size(good_dots_stats(1).Centroid,2) == 3 %3D image
                % Readjust dots centroid to reference frame before autocrop
                for l=1:dots(k).counts
                    dots(k).properties(l).Centroid = dots(k).properties(l).Centroid + [x_ind(1)-1 y_ind(1)-1 0 ];
                end
                centroids{k} = reshape([dots(k).properties.Centroid], 3, length(dots(k).properties))';
            else
                
                % Readjust dots centroid to reference frame before autocrop
                for l=1:dots(k).counts
                    dots(k).properties(l).Centroid = dots(k).properties(l).Centroid + [x_ind(1)-1 y_ind(1)-1 ];
                end
                
                centroids{k} = reshape([dots(k).properties.Centroid], 2, length(dots(k).properties))';
                centroids{k}(:,3)=1; % on plane 1 in 3D
            end
        end
    end
    
    % compare other channels to first channel and eliminate dots
    % that colocalize with dots in the first channel
    OVERLAP_SIZE = 2; channels_to_compare = 2:num_channels; channel_standard = 2;  
    dots_std_to_remove = []; 
%     for k = channels_to_compare(channels_to_compare ~= channel_standard)
%         if ~isempty(centroids{channel_standard}) && ~isempty(centroids{k})
%             [pairwise_distances, indices] = pdist2(centroids{2}, centroids{k}, 'euclidean', 'Smallest', 1);
%             dots_to_keep = find(~ismember(1:length(dots(k).properties), find(pairwise_distances < OVERLAP_SIZE)));
%             dots(k).properties = dots(k).properties(dots_to_keep); dots(k).counts = length(dots_to_keep);
%             dots_std_to_remove = [unique(indices(pairwise_distances < OVERLAP_SIZE))'; dots_std_to_remove(:)];
%         end
%     end
%     
    dots_std_to_keep = find(~ismember(1:length(dots(channel_standard).properties), dots_std_to_remove));
    dots(channel_standard).properties = dots(channel_standard).properties(dots_std_to_keep); dots(channel_standard).counts = length(dots_std_to_keep);
end


% --------------------------------------------------------------------
function lapFrame = logMask(im)   %**check that have the right logMask function

    k = [-4 -1  0 -1 -4;...
         -1  2  3  2 -1;...
          0  3  4  3  0;...
         -1  2  3  2 -1;...
         -4 -1  0 -1 -4];

    lapFrame = imfilter(im,k,'repl');
end

function [cropped_image,y_ind,x_ind]= autocrop(image)
        sum_image=sum(image,3);
            x_min = find(sum(sum_image,1)>0,1,'first');
            x_max = find(sum(sum_image,1)>0,1,'last');
            
            y_min = find(sum(sum_image,2)>0,1,'first');
            y_max = find(sum(sum_image,2)>0,1,'last');
        
        y_ind = max(1,y_min-10):min(size(image,1),y_max+10);
        x_ind = max(1,x_min-10):min(size(image,2),x_max+10);
        cropped_image  = image(y_ind,x_ind,:);
end