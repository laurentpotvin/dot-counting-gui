function convert_stack_to_csv(stack)
[filepath,name,ext] = fileparts(stack);
dotfile = strcat(filepath,'/',name,'_dots.csv');
areafile = strcat(filepath,'/',name,'_area.csv');
data = load(stack);
[dots, dot_sizes, frame_info, cell_area] = extract_dot_info(data);
csvwrite(dotfile, dots);
csvwrite(areafile, cell_area);
end