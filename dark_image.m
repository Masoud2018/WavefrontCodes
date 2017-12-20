clear all; close all;

this_file = [mfilename '.m'];
%matlab_path_1 = '~/Documents/software/matlab/';
matlab_path_2 = '~/FLASH_2015/DATA';

%addpath([matlab_path_1 'external_routines/export_fig/']);
addpath(matlab_path_2);
addpath([matlab_path_2 '/function_library']);

text_font_size = 18;

set(0, ...
    'defaultaxesfontsize',   text_font_size, ...
    'defaultaxeslinewidth',   2, ...
    'defaultlinelinewidth',   2.5, ...
    'defaultpatchlinewidth',  2, ...
    'defaultaxesfontname', 'Arial', ...
    'defaultaxesbox', 'on');

% avoid changing black background to white upon printing
set(0,'defaultfigureInvertHardCopy','Off');

% plot initialization command
cmd_init_plot = ['fig_ind = fig_ind+1; figure(fig_ind);',...
'set(gcf,''color'',''white'',''units'',''normalized'',''position'',[0.1 0.1 0.6 0.8]);'];

fig_ind = 0;

% new colormap darkjet (Robin's live viewer)
darkjet = jet(256);
%darkjet(1,:) = 0;

invgray = gray(256);
invgray = invgray(256:-1:1,:);


%% dark preparation 

% start time
time1 = '2015-06-19 06:56:21';

% end time 
time2 = '2015-06-19 07:01:09';

% issue system command

[statusc,resultc] = system(['/Users/mehrjoom/FLAsH_2015/DATA/function_library/which_images_in_time.sh ' ...
    '"' time1 '" "' time2 '" | sort > /Users/mehrjoom/FLASH_2015/DATA/analysis/dark_images/file_list_dark_19.txt']);

% save dir
save_dir = '/Users/mehrjoom/FLASH_2015/DATA/analysis/dark_images';
 
% save analysis?
  save_analysis = true;

 % get file list
 fid = fopen('/Users/mehrjoom/FLASH_2015/DATA/analysis/dark_images/file_list_dark_19.txt');
 files = textscan(fid,'%s');
 fclose(fid);
 files = files{1};
 P.N_files = length(files);
 
 % read test frame
 test = read_raw(files{1});
 dark = zeros(size(test,1),size(test,2));
 dark_1d = zeros(P.N_files,1);

%% read data

% filename of image to read
eval(cmd_init_plot);
for n = 1:P.N_files
    dark_single = read_raw(files{n});
    dark = dark + dark_single;
    dark_1d(n) = sum(dark_single(:));
    
    if ~mod(n,1)
        title_str = sprintf(['Image: ' files{n}(end-32:end)]);
        imagesc(rot90(dark_single)); axis equal tight; colorbar; colormap(jet); title(title_str,'Interpreter','none');
        zoom(5); drawnow;
    end
end

dark = rot90(dark)/P.N_files;
dark_19 = double(dark) ;
%% parameters

% threshold for noise cutoff
P.noise_threshold = 0;

% name of data file
P.fname = 'ANDOR.CAM1_20150622T172402_73_779918721.raw';

% name of dark file
P.dname = 'ANDOR.CAM1_20150622T171253_383_779912036.raw';

% save dir
P.save_dir = '/Users/mehrjoom/FLASH_2015/DATA/analysis/';

% data dir
P.data_dir = '/Users/mehrjoom/local_folder2/data_FLASH2015/raw_data/ANDOR_RAW/';

% save analysis?
P.save_analysis = true ;

%source to aperture distance (m)
P.z = 71.15 ;
% distance aperture-mirror center (m)
P.z01 = 3.85 ;
%P.z01 = 17 ;
% focal length (m) 
P.z12 = 2 ;
% distance focus-to-detector (m)
P.z23 = 1.2;
% photon wavelength (m)
P.lambda = 15*1e-9;
% wave number (m^(-1))
P.k = 2 * pi / P.lambda ;

% method for finding beam center ('cm' or 'manual')
P.method_beam_center = 'cm';

% size of the image for reconstruction
P.Ny = 2048;
P.Nx = 2048;

% Detector plane (z3)
% pixel width (usually detector pixel width)
P.dx3 = 13e-6;
P.dy3 = 13e-6;
% 2D grid
[P.X3,P.Y3] = meshgrid(((1:P.Nx)-floor(P.Nx/2)-1)*P.dx3,((1:P.Ny)-floor(P.Ny/2)-1)*P.dy3);
X3_axis = ((1:P.Nx)-floor(P.Nx/2)-1)*P.dx3;
Y3_axis = ((1:P.Ny)-floor(P.Ny/2)-1)*P.dy3;

% Focal plane (z2)
% pixel widths
P.dx2 = P.lambda*P.z23/(P.Nx*P.dx3);
P.dy2 = P.lambda*P.z23/(P.Ny*P.dy3);
% 2D grid
[P.X2,P.Y2] = meshgrid(((1:P.Nx)-floor(P.Nx/2)-1)*P.dx2,((1:P.Ny)-floor(P.Ny/2)-1)*P.dy2);
X2_axis = ((1:P.Nx)-floor(P.Nx/2)-1)*P.dx2;
Y2_axis = ((1:P.Ny)-floor(P.Ny/2)-1)*P.dy2;

% Pupil plane (z1) / mirror entrance plane
% pixel widths
P.dx1 = P.lambda*P.z12/(P.Nx*P.dx2);
P.dy1 = P.lambda*P.z12/(P.Ny*P.dy2);
% 2D grid
[P.X1,P.Y1] = meshgrid(((1:P.Nx)-floor(P.Nx/2)-1)*P.dx1,((1:P.Ny)-floor(P.Ny/2)-1)*P.dy1);
X1_axis = ((1:P.Nx)-floor(P.Nx/2)-1)*P.dx1;
Y1_axis = ((1:P.Ny)-floor(P.Ny/2)-1)*P.dy1;

% circular entrance pupil
P.Dp = 25e-3;

% Circular Aperture
P.Da = 10e-3 ;


fprintf('dx1 = %3.1f um.\n',P.dx1*1e6);
fprintf('dx2 = %3.1f nm.\n',P.dx2*1e9);
fprintf('dx3 = %3.1f um.\n',P.dx3*1e6);
fprintf('Diffraction limited focal spot diameter (FWHM) for aperture diameter %3.3f mm is %3.3f um\n.',P.Dp*1e3,1.03*P.lambda*P.z12/P.Dp*1e6);

%% read data

title_str = sprintf(['Image: ' P.fname  '\n Dark: ' P.dname]);

% filename of dark image to read
if numel(P.dname) > 0
    dark = read_raw([P.data_dir P.dname]);
else
    dark = 0;
end
dark = (rot90(dark));
%save('/Users/mehrjoom/FLASH_2015/DATA/analysis/DARK.mat','dark');

%filename of image to read
data = read_raw([P.data_dir P.fname]);
data = (rot90(data));

data = max(data - dark,P.noise_threshold);

%%

% make valid mask
P.I_valid = ones(size(data));
P.I_valid(1:4,:) = 0;

switch P.method_beam_center

    case 'cm'
   
        [P.cx,P.cy] = center_of_mass(data.*P.I_valid);        
        CX = round(P.cx); CY = round(P.cy);

    case 'manual'

        % take the values that have been defined in the parameter section
        
end

% embed the data

P.I_exp = embed_iamge(zeros(P.Ny,P.Nx),data.*P.I_valid,CX,CY);
P.I_valid = embed_iamge(zeros(P.Ny,P.Nx),P.I_valid,CX,CY);

eval(cmd_init_plot);
imagesc(X3_axis,Y3_axis,log10(data)); axis equal tight; colorbar; colormap(invgray);
title(title_str,'Fontsize',12,'Interpreter','none');
im_info = ['Max.: ' num2str(max(data(:)),'%2.2e') ', Min.: ' num2str(min(data(:)),'%2.2e')];
im_info = [im_info ', Sum: ' num2str(sum(data(:)),'%2.2e')];
text(0.02,0.02,im_info,'units','normalized','Fontsize',12,'Color','k');

eval(cmd_init_plot);
imagesc(X3_axis,Y3_axis,P.I_exp); axis equal tight; colorbar; colormap(invgray);
title(title_str,'Fontsize',12,'Interpreter','none');
im_info = ['Max.: ' num2str(max(P.I_exp(:)),'%2.2e') ', Min.: ' num2str(min(data(:)),'%2.2e')];
im_info = [im_info ', Sum: ' num2str(sum(P.I_exp(:)),'%2.2e')];
text(0.02,0.02,im_info,'units','normalized','Fontsize',12,'Color','k');

eval(cmd_init_plot);
imagesc(P.I_valid); axis equal tight; colorbar; colormap(invgray);
title(title_str,'Fontsize',12,'Interpreter','none');

%% reconstruction

P.sigma = 40;
P.alpha = 1;
P.z_source = 50;

P.Amp_exp = sqrt(P.I_exp);
P.psi = zeros(P.Ny,P.Nx);

G = fspecial('gaussian', 5*P.sigma, P.sigma);
phase_0 = imfilter(rand(P.Ny,P.Nx),G,'circular');
phase_0 = (phase_0-min(phase_0(:)))/(max(phase_0(:)-min(phase_0(:))));


% initial guess (random phases)
P.Q = P.Amp_exp.*exp(1i*P.alpha*2*pi*phase_0);
P.Q = P.Q.*exp(1i*pi*(P.X3.^2+P.Y3.^2)/(P.lambda*P.z_source));


%%
% % initialize pupil function
% P.P = double(P.support);
% % normalize to experimental photon count
% P.P = P.P * norm(P.Amp_exp(:).*P.I_valid(:))/norm(P.P(:));

% number of iterations
P.N_it = 10;

P.err = zeros(P.N_it,1);
P.err_denom = sum(P.Amp_exp(:).*P.I_valid(:)).^2;
for ii=1:P.N_it
    tic;
    P.it = ii;
    P = iterate_wavefield(P);
    P.err(ii) = P.err_it;
    fprintf('Fourier space error: %3.3f.\n',P.err_it);
    fprintf('Iteration %i of %i.\n',ii,P.N_it);
    toc;
end

%% plot result

% %pupil function
pupil = prop_nf_pa(P.P,P.lambda,-P.z01,P.dx1,P.dy1);

eval(cmd_init_plot);
imagesc(hsv2rgb(im2hsv(pupil*exp(1i*1)/max(abs(pupil(:))),1))); axis equal tight;

eval(cmd_init_plot);
imagesc(X2_axis,Y2_axis,hsv2rgb(im2hsv(P.psi*exp(1i*1)/max(abs(P.psi(:))),1))); axis equal tight;

eval(cmd_init_plot);
imagesc(hsv2rgb(im2hsv(P.Q*exp(1i*1)/max(abs(P.Q(:))),1))); axis equal tight;

eval(cmd_init_plot);
imagesc(abs(P.Q/max(abs(P.Q(:))))); axis equal tight;

%eval(cmd_init_plot);
%imagesc(abs(P.Q/max(abs(P.Q(:))))); axis equal tight;

%% back propagation of the illumination function

N = 60;
dz = 2e-3;

z_range = ((1:N)-floor(N/2)-1)*dz;
PROBE = zeros(P.Ny,P.Nx,N);
sharpness = zeros(N,1);

for jj = 1:N
   PROBE(:,:,jj) = prop_nf_pa(P.psi,P.lambda,z_range(jj),P.dx2,P.dy2,'quiet');
   sharpness(jj) = sum(sum(abs(squeeze(PROBE(:,:,jj))).^4));   
   
   %imagesc(hsv2rgb(im2hsv(PROBE(:,:,jj)/max(max(abs(PROBE(:,:,jj)))),1))); 
   imagescg(abs(PROBE(:,:,jj)/max(max(abs(PROBE(:,:,jj))))),1); 
   axis tight; set(gca,'visible','off');
   axlim = axis(gca);
   text(0.7*(axlim(2)-axlim(1))+axlim(1),0.95*(axlim(4)-axlim(3))+axlim(3),['z = ' num2str(z_range(jj)*1e3,'%4.3f') ' mm'],'FontSize',text_font_size,'Color','w');
   %drawnow;

   if ~mod(jj,10)
       fprintf('%2.2f per cent finished.\n ',jj/N*100);
   end

end


%% Plot through-focus-propagation

% coordinate axes (real space), for imagesc
% end points
X_1 = ([1,P.Nx]-floor(P.Nx/2)-1)*P.dx1;
Y_1 = ([1,P.Ny]-floor(P.Ny/2)-1)*P.dy1;
% full axis vectors
X_1_axis = ((1:P.Nx)-floor(P.Nx/2)-1)*P.dx1;
Y_1_axis = ((1:P.Ny)-floor(P.Ny/2)-1)*P.dy1;


% Focal plane estimation

eval(cmd_init_plot);
plot(z_range*1e3,sharpness/max(sharpness));
xlabel('z [mm]'); ylabel('Sharpness [a.u.]');

z_focus = z_range(sharpness == max(sharpness));

% plane to be used for 2D complex plot of the probe
z_select = z_focus;
% manual focus selection
%z_select = z_focus;

% zy-plane
eval(cmd_init_plot); set(gcf,'color','black'); % make bg black
% normalized
max_PROBE = max(max(abs(squeeze(PROBE(:,floor(P.Nx/2)+1,:))).^2));
imagesc([z_range(1),z_range(end)]*1e3,Y_1*1e6,abs(squeeze(PROBE(:,floor(P.Nx/2)+1,:))).^2/max_PROBE);
%load('klaus_colormaps','logmap'); set(fig_ind,'Colormap',logmap); 
cb = colorbar; set(get(cb,'ylabel'),'String','Intensity [norm. units]','Fontsize',text_font_size);
set(gca,'XColor','w','YColor','w','XMinorTick','On','YMinorTick','On');
xlabel('z [mm]','Color','w'); ylabel('y [\mum]','Color','w');
% plot sample plane/reconstruction plane
line([0 0], Y_1*1e6 , 'color','white','linestyle',':');
% plot determined focal plane
line([z_select*1e3 z_select*1e3], Y_1*1e6 , 'color','white','linestyle','--');
l = legend({'Sample','Focus'}); set(l,'Color','none','TextColor','w');
% plot pinhole plane
%line([z_pinhole*1e3 z_pinhole*1e3], Y_1*1e6 , 'color','white','linestyle','-.');

% zx-plane
eval(cmd_init_plot); set(gcf,'color','black'); % make bg black
% normalized
max_PROBE = max(max(abs(squeeze(PROBE(floor(P.Ny/2)+1,:,:))).^2));
imagesc([z_range(1),z_range(end)]*1e3,X_1*1e6,abs(squeeze(PROBE(floor(P.Ny/2)+1,:,:))).^2/max_PROBE);
%load('klaus_colormaps','logmap'); set(fig_ind,'Colormap',logmap); 
cb = colorbar; set(get(cb,'ylabel'),'String','Intensity [norm. units]','Fontsize',text_font_size);
set(gca,'XColor','w','YColor','w','XMinorTick','On','YMinorTick','On');
xlabel('z [mm]'); ylabel('x [\mum]');
% plot sample plane/reconstruction plane
line([0 0], X_1*1e6 , 'color','white','linestyle',':');
% plot determined focal plane
line([z_select*1e3 z_select*1e3], X_1*1e6 , 'color','white','linestyle','--');
l = legend({'Sample','Focus'}); set(l,'Color','none','TextColor','w');
% plot pinhole plane
%line([z_pinhole*1e3 z_pinhole*1e3], X_1*1e6 , 'color','white','linestyle','-.');

%% propagation to focus

eval(cmd_init_plot);
% reset geometry a bit
set(gcf,'color','white','units','normalized','position',[0.1 0.1 0.7 0.8]);
Probe_select = prop_nf_pa(P.psi,P.lambda,z_select,P.dx2,P.dy2,'quiet');
imagesc(X2_axis,Y2_axis,hsv2rgb(im2hsv(Probe_select/max(abs(Probe_select(:))),1))); 
axis equal tight; set(gca,'visible','off');
add_scalebar(gca,0.5e-6,P.Nx*P.dx2,'nm','Fontsize',text_font_size);
axlim = axis(gca);
text(0.7*(axlim(2)-axlim(1))+axlim(1),0.95*(axlim(4)-axlim(3))+axlim(3),['z = ' num2str(z_focus*1e3,'%4.3f') ' mm'],'FontSize',text_font_size,'Color','w');
%add_colorwheel(gca,0);

 fig_ind = fig_ind+1; figure(fig_ind); set(gcf,'color','white','units','normalized','position',[0.1 0.1 0.4 0.6]);
 imagesc(hsv2rgb(im2hsv(Probe_select*exp(1i*1)/max(abs(Probe_select(:))),1))); axis equal tight;
 add_scalebar(gca, 500e-9, P.Nx*P.dx2, 'nm');
title('Focal plane');

int_Probe = abs(Probe_select).^2/max(abs(Probe_select(:)).^2);
% subscript indices of the global maximum
[ind_max_i,ind_max_j] = find(int_Probe == max(int_Probe(:)));

fig_ind = fig_ind+1; figure(fig_ind); set(gcf,'color','white','units','normalized','position',[0.1 0.1 0.4 0.6]);
% reset geometry a bit
set(gcf,'color','white','units','normalized','position',[0.1 0.1 0.6 0.8]);
% determine the FHWM in horizontal direction
%[FWHM_x,x1,x2,HM] = get_FWHM_1D(X2_axis*1e6,int_Probe(ind_max_i,:));
plot(X2_axis*1e6,int_Probe(ind_max_i,:)); %xlim(X_1*1e6);
%[figx figy] = dsxy2figxy(gca,[x1 x2],[HM HM]);
%annotation(gcf,'DoubleArrow',figx,figy);
xlabel('x [\mum]'); ylabel('Intensity [norm. units]');
%text(x2 + 2.5*x2,HM,['FWHM = ' num2str(FWHM_x*1e3,'%3.0f') ' nm'],'FontSize',text_font_size,'Rotation',90,'HorizontalAlignment','Center');

fig_ind = fig_ind+1; figure(fig_ind); set(gcf,'color','white','units','normalized','position',[0.1 0.1 0.4 0.6]);
% reset geometry a bit
set(gcf,'color','white','units','normalized','position',[0.1 0.1 0.6 0.8]);
% determine the FHWM in vertical direction
%[FWHM_y,y1,y2,HM] = get_FWHM_1D(Y2_axis*1e6,int_Probe(:,ind_max_j));
plot(Y2_axis*1e6,int_Probe(:,ind_max_j)); %xlim(Y_1*1e6);
%[figx figy] = dsxy2figxy(gca,[y1 y2],[HM HM]);
%annotation(gcf,'DoubleArrow',figx,figy);
xlabel('y [\mum]'); ylabel('Intensity [norm. units]');
%text(x2 + 2.5*x2,HM,['FWHM = ' num2str(FWHM_y*1e3,'%3.0f') ' nm'],'FontSize',text_font_size,'Rotation',90,'HorizontalAlignment','Center');

%fprintf('Intensity FWHM_x = %3.0f nm\n',FWHM_x*1e3);
%fprintf('Intensity FWHM_y = %3.0f nm\n',FWHM_y*1e3);
%fprintf('FWHM_y/FWHM_x = %2.2f\n',FWHM_y/FWHM_x);






