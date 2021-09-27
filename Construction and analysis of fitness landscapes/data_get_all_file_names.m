function file_list = data_get_all_file_names(folder)
% get the names of all data files in a folder
if nargin < 1
    folder = 'reads_parsed'
end

d = dir(folder);
file_list = {d([3:end]).name};
for fi = 1:numel(file_list)
    file_list{fi} = sprintf('%s/%s', folder, file_list{fi}); % add on path to the folder
end
    
    