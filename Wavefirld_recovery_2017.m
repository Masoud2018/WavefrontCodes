
% This code retrieves the wave fields measured at the FLASH facility 
% at the beamline BL2 for the specific wavelength and geometry. 
% For further applications, the iterative algorithm may change thoroughly
% depending on the geometry and the wavelength. The criteria for the
% sampling have been described in Masoud's thesis on pages 52-53.  


%% Start

close all, clear all;

%% Paths (sunject to change)

this_file = [mfilename '.m'];
matlab_path_3 = 'H:\phd_first_step\jobs\main_FCDI_codes\codes' ;
matlab_path_4 = 'H:\phd_first_step\jobs\main_FCDI_codes' ;
addpath(matlab_path_3);
addpath(matlab_path_4);
addpath([matlab_path_3 '\function_library_2']);
path = [matlab_path_3 '\function_library_2'] ;

%% Important Data 

% Here specify the parameter follow:
% number of iterations
it = 20 ;
% Pixel Number
N = 1024 ;
% Wave lenght\number
lambda = 14.7*1e-9 ;
k = 2 * pi / P.lambda ;
% Axis's upscaling factor
R = 1e3 ;
RR = 1e6 ;
% shift intensity's mask threshold
limit = 0.11 ;
threshold = 0.15 ;
% Reconstruction parameters
N_start = 1 ;
N_end = 200/2 ;
gap = 1 ;
% initialization parameter
sigma = 5;
alpha = 1 ;
% erorr measure
tmr = 1e6 ;
% method for finding beam center ('cm' or 'manual')
method_beam_center = 'cm';


%% Set the parameters of the pulses ('A' '1-4' '10-5-3') 

stamp_0 = num2str(3) ;
stamp_1 = num2str(5) ;
pulses = ['A-' stamp_0 '-' stamp_1 ] ;

[Regime,regime,Undulator,undulator,RAP,O,AP,APS] = pulse(pulses) ;

%% Specify the geomerty 

% detector plane(m)
dx4 = 13e-6;
dy4 = 13e-6;
% 2D grid
[X4,Y4] = meshgrid(((1:N)-floor(N/2)-1)*dx4,((1:N)-floor(N/2)-1)*dy4);
X4_axis = ((1:N)-floor(N/2)-1)*dx4;
Y4_axis = ((1:N)-floor(N/2)-1)*dy4;

% Lens_Detector Setting
Distance = 3.2 ;
z23 = 2.054 ;
z34 = Distance - z23 ;

% image plane
dx3 =  lambda*z34/(N*dx4) ;
dy3 =  lambda*z34/(N*dx4) ;
% 2D grid
[X3,Y3] = meshgrid(((1:N)-floor(N/2)-1)*dx3,((1:N)-floor(N/2)-1)*P.dy3);
X3_axis = ((1:N)-floor(N/2)-1)*dx3;
Y3_axis = ((1:N)-floor(N/2)-1)*dy3;


% lens plane
dx2 = lambda*z23/(N*dx3);
dy2 = lambda*z23/(N*dy3);
% 2D grid
[X2,Y2] = meshgrid(((1:N)-floor(N/2)-1)*dx2,((1:N)-floor(N/2)-1)*dy2);
X2_axis = ((1:N)-floor(N/2)-1)*dx2;
Y2_axis = ((1:N)-floor(N/2)-1)*dy2;



%  aperture
% pixel widths
dx1 = dx2 ; 
dy1 = dy2 ;
% 2D grid
[X1,Y1] = meshgrid(((1:N)-floor(N/2)-1)*dx1,((1:N)-floor(N/2)-1)*dy1);
X1_axis = ((1:N)-floor(N/2)-1)*dx1;
Y1_axis = ((1:N)-floor(N/2)-1)*dy1;



% aperture-source-lens
z01 = 71.15 ;
f = 2 ;
z12 = 3.85 ;
ap = circ(X1,Y1,O) ;
%% Fourier Domain Configuration

du = 2*pi / (N * dx1) ;
[ux,uy] = meshgrid(du *(-N/2:N/2-1),du * (-N/2:N/2-1)) ;
Kernel_Factor = exp(-1i*(ux.^2+uy.^2)*z12/(2*k));
Invers_Kernel_Factor = conj(Kernel_Factor) ;

%% dark preparation (Subject to cahnge)

load([matlab_path_4 '\Result\Raw_Text_Dark\dark_A_' stamp_0 '_' APS '.mat']) ;

%% data preparation (Subject to change)
 
 % get file list
 fid = fopen([matlab_path_4 ...
    '\Result\Raw_Text_Shot\file' '_' APS '_' regime '_' undulator ...
    '.txt']);
 files = textscan(fid,'%s');
 fclose(fid);
 files = files{1};
 P.N_files = length(files);
 
%% read data

% filename of image to read
Data_sum = zeros(N,N) ;
for n = 1 : N_end
    Data(:,:,n) = read_raw(files{n});
    data = squeeze(rot90(Data(:,:,n)))  ;
    imagesc(data),axis equal tight,title(num2str(n)),drawnow

end

%% consistent range

N_se = 15 ;
N_ne = 15 ;
Data_best(:,:,1:(N_ne-N_se)+1) = Data(:,:,N_se:N_ne) ;
count = N_ne-N_se ;

%% Reconstruction

tic
size = 10e-3 ;
for ii = 1 : gap : (N_ne-N_se)+1

    %%%% data preparation %%%%%%%%%   
  index = num2str(ii) ;
  data = squeeze(rot90(Data_best(:,:,ii)))  ;    

  %%%%% Save Pured Diffraction Pattern %%%%%%
  
% back ground subtraction
   data = max(data - dark,threshold);
% centering the image
   data = shift_intensity(data,limit,N) ;


%%%%%%% Initialization%%%%%%
  Psi_initial = FCDI_intialization(data,sigma,alpha,N) ;
  I_exp = data / max(data(:)) ;
  E_exp = sqrt(I_exp) ;
  
%%%%%% Iterative algorithm %%%%%%
  [E_aperture,E_focus,E_det,Error,Support] = iterative_function(Psi_initial,Invers_Kernel_Factor,Kernel_Factor ...
    ,E_exp,it,lambda,k,z34,z23,X2,Y2,X3,Y3,f) ;

fprintf('The last iteration:%2d\n',ii) ;
end

time = toc/60 ;
fprintf('Elapsed time to analyze each images:%3.1f',time/count) ;



