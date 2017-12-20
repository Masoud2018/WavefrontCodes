function f = count_image(time1,time2)

matlab_path_3 = '~/FLASH_FINAL_ANALYSIS/codes' ;
matlab_path_4 = '~/FLASH_FINAL_ANALYSIS/result/Raw_Text' ;
% issue system command

[status,result] = system([matlab_path_3 ...
    '/function_library_2/which_images_in_time.sh ' ...
    '"' time1 '" "' time2 '" | sort >' matlab_path_4 ...
    '/images.txt']);

 
 % get file list
 fid = fopen([matlab_path_4 ...
    '/images.txt']);
 files = textscan(fid,'%s');
 fclose(fid);
 files = files{1};
 f = length(files);
 
end