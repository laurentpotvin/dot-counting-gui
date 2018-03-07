function [dots, dot_sizes, frame_info] = extract_dot_info(data)
dots = []; %zeros(size(data.stack.frame, 2), length([data.stack.frame(1).cell(:).dots.counts]));
frame_info = [];
number_channels = -1;
% Find the number of channels
for k = 1:size(data.stack.frame, 2)
    frame_data = data.stack.frame(k);
    for p = 1:length(frame_data.cell)
        if isfield(frame_data.cell(p), 'dots')
            if length(frame_data.cell(p).dots) > number_channels
                number_channels = length(frame_data.cell(p).dots);
            end
        end
    end
end

if number_channels  >0
    dot_sizes = cell(1,number_channels);
    
    
    for k = 1:size(data.stack.frame, 2)
        frame_data = data.stack.frame(k);
        for p = 1:length(frame_data.cell)
            if isfield(frame_data.cell(p), 'dots')
                frame_info(end + 1) = k;
                dots(end + 1, :) = [frame_data.cell(p).dots.counts];
                for m = 1:length(frame_data.cell(p).dots)
                    temp = dot_sizes{m};
                    temp = cat(1, temp, frame_data.cell(p).dots(m).properties.Area);
                    dot_sizes{m} = temp;
                end
            end
        end
    end
end

end