%% dark preparation 

% start time
time1 = '2015-06-22 17:12:51';

% end time 
time2 = '2015-06-22 17:16:53';

% issue system command

[status_dark,result_dark] = system([matlab_path_3 ...
    '/function_library_2/which_images_in_time.sh ' ...
    '"' time1 '" "' time2 '" | sort >' matlab_path_4 ...
    '/result/Raw_Text/file' '_' APS '_' regime '_' undulator '.txt']);

 % get file list
 fid = fopen([matlab_path_4 ...
     '/result/Raw_Text/file' '_' APS '_' regime '_' undulator '.txt']);
 files = textscan(fid,'%s');
 fclose(fid);
 files = files{1};
 P.N_files = length(files);
 
 % read test frame
 test = read_raw(files{1});
 dark = zeros(size(test,1),size(test,2));
 dark_1d = zeros(P.N_files,1);

%% read dark data

% filename of image to read
for n = 1:P.N_files
    dark_single = read_raw(files{n});
    dark = dark + dark_single; 
end

dark = rot90(dark)/P.N_files;
dark = double(dark) ;
dark_n = dark(:) ;
mean_dark = mean(dark_n) ;
std_dark = std2(dark) ;
%%
fprintf('mean of dark:%4.4f\t\n',mean_dark);
fprintf('Standard Deviation of dark:%4.4f\t\n',std_dark);
fprintf('std/mean:%4.4f\t\n',std_dark/mean_dark);








