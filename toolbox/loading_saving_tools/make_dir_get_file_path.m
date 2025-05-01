function [file_path] = make_dir_get_file_path(dir_main, dir_name, struct_file)

 dir = fullfile(dir_main, dir_name);
 if ~exist(dir, 'dir')
     mkdir(dir);
 end

 file_path = fullfile(dir,struct_file)