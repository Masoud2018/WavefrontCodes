

close all, clear all;

%% Paths

this_file = [mfilename '.m'];
matlab_path_3 = 'H:\phd_first_step\jobs\main_FCDI_codes\codes' ;
matlab_path_4 = 'H:\phd_first_step\jobs\main_FCDI_codes' ;
matlab_path_5 = 'D:\partial_coherence_study' ;

%%

% Set the parameters of the pulses ('A' '1-4' '10-5-3')
stamp_0 = num2str(2) ;
stamp_1 = num2str(3) ;
pulses = ['A-' stamp_0 '-' stamp_1 ] ;

[Regime,regime,Undulator,undulator,RAP,O,AP,APS] = pulse(pulses) ;

[time1,time2] = dark_frames(Undulator,Regime) ;
[time3,time4] = data_frames(Undulator,Regime,AP) ;

%%
%%% Unix path
matlabpath_0 = 'H:/phd_first_step/jobs/main_FCDI_codes/codes' ;
matlabpath_1 = 'H:/phd_first_step/jobs/main_FCDI_codes' ;

[status_dark,result_dark] = system(['C:\cygwin64\bin\bash --login ' matlabpath_0 ...
    '/function_library_2/Time_Frame.sh' ' "' time1 '" "' time2 '" | sort >'...
    matlabpath_1 '/Result/Raw_Text_Dark/file' '_' APS '_' regime '_' undulator '.txt']);

 % get file list
 fid = fopen([matlab_path_4 ...
     '\Result\Raw_Text_Dark\file' '_' APS '_' regime '_' undulator '.txt']);
 files = textscan(fid,'%s');
 fclose(fid);
 files = files{1};
 P.N_files = length(files);
 
 % read test frame
 test = read_raw(files{1});
 dark = zeros(size(test,1),size(test,2));
 dark_1d = zeros(P.N_files,1);

%% read dark data and average them

% filename of image to read
for n = 1:P.N_files
    dark_single = read_raw(files{n});
    dark = dark + dark_single; 
    %dark_1d(n) = sum(dark_single(:));
end

dark = rot90(dark)/P.N_files;
dark = double(dark) ;

save([matlab_path_4 '\Result\Raw_Text_Dark\dark_A_' stamp_0 '_' APS ,'.mat'],'dark') ;
%% data preparation 


% issue system command

[status_data,result_data] = system(['C:\cygwin64\bin\bash --login ' matlabpath_0 ...
    '/function_library_2/Time_Frame.sh' ' "' time3 '" "' time4 '" | sort >'...
    matlabpath_1 '/Result/Raw_Text_Shot/file' '_' APS '_' regime '_' undulator '.txt']);







