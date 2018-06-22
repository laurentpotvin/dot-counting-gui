function change_stack_basepath(stack_name, old_path, new_path)
% function change_stack_basepath(stack_name, old_path, new_path)
% Change the basepath for a stack
% stack_name: filename for loading the stack
% old_path: basepath that needs to be replaced (e.g. /path/to/old/files/)
% new_path: basepath where the files are currently located (e.g. /path/to/new/files/)
% Note: make sure that both old_path and new_path ends with a slash (or
% not)
% The new stack is saved with a new_path extension
load(stack_name);

old_path_length = strlength(old_path);

if strncmp(old_path, stack.name, old_path_length)
    stack.name = [new_path stack.name(old_path_length+1:end)];
else
    error('cannot find old path')
end
for i=1:length(stack.image_path_cell)
    if iscell(stack.image_path_cell{i}) %multiple hyb
       for j=1:length( stack.image_path_cell{i})
           stack.image_path_cell{i}{j} = [new_path stack.image_path_cell{i}{j}(old_path_length+1:end) ];
           
       end
    else
        stack.image_path_cell{i} = [new_path stack.image_path_cell{i}(old_path_length+1:end) ];
    end
end

save([stack_name 'new_path'],'stack')
end