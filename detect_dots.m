function dots = detect_dots(imdata, seg_im, num_channels, thresholds)
    
    dots = struct; DOT_AREA = 2; centroids = cell(num_channels, 1);
    for k = 1:num_channels
        im_stack = imdata{k};
        for p = 1:size(im_stack, 3)
            bg_sub_imstack(:, :, p) = im_stack(:, :, p) - medfilt2(im_stack(:, :, p), [5 5]); %logMask(im_stack(:, :, p)).*uint16(logical(seg_im));
        end
        log_filt_imstack = bg_sub_imstack; %imgaussfilt3(bg_sub_imstack, 0.8); %0*im_stack;
        pot_dots_im = imregionalmax(log_filt_imstack);
                
        thresh_stack = bsxfun(@times, bg_sub_imstack, uint16(seg_im)); 
        thresh_stack(thresh_stack < thresholds(k)) = 0;
        
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
                centroids{k} = reshape([good_dots_stats.Centroid], 3, length(good_dots_stats))';
            else
                centroids{k} = reshape([good_dots_stats.Centroid], 2, length(good_dots_stats))';
                centroids{k}(:,3)=1; % on plane 1 in 3D
            end
        end
    end
    
    % compare other channels to first channel and eliminate dots
    % that colocalize with dots in the first channel
    OVERLAP_SIZE = 2; channels_to_compare = 2:num_channels; channel_standard = 2;  
    dots_std_to_remove = []; 
    for k = channels_to_compare(channels_to_compare ~= channel_standard)
        if ~isempty(centroids{channel_standard}) && ~isempty(centroids{k})
            [pairwise_distances, indices] = pdist2(centroids{2}, centroids{k}, 'euclidean', 'Smallest', 1);
            dots_to_keep = find(~ismember(1:length(dots(k).properties), find(pairwise_distances < OVERLAP_SIZE)));
            dots(k).properties = dots(k).properties(dots_to_keep); dots(k).counts = length(dots_to_keep);
            dots_std_to_remove = [unique(indices(pairwise_distances < OVERLAP_SIZE))'; dots_std_to_remove(:)];
        end
    end
    
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