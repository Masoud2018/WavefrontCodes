close all, clear all,
matlab_path_2 = '~/FLASH_2015/DATA';
% object to aperture  (m) 
P.z12 = 71.15 ;
% Aperture to mirror (m)
P.z23 = 3.85 ;
% Mirror to detector
MTD = 3.234 ;
% mirror to image  (m)
%P.z34 = 2.0548 ;%Theoretically the image position
P.z34 = 2.0;% an emprically measured position(can be foci or around it)
% image to detector plane  (m)
P.z45 = MTD - P.z34 ;
% focal length
P.f = 2.0 ;
% photon wavelength (m)
P.lambda = 15e-9;
% wave number (m^(-1))
P.k = 2 * pi / P.lambda ;
% size of the image for reconstruction
P.Ny = 256;
P.Nx = 256;

% aperture size(m)
P.O = 2e-3;

% detector plane(m)
P.dx4 = 13e-6;
P.dy4 = 13e-6;
% 2D grid
[P.X4,P.Y4] = meshgrid(((1:P.Nx)-floor(P.Nx/2)-1)*P.dx4,((1:P.Ny)-floor(P.Ny/2)-1)*P.dy4);
X4_axis = ((1:P.Nx)-floor(P.Nx/2)-1)*P.dx4;
Y4_axis = ((1:P.Ny)-floor(P.Ny/2)-1)*P.dy4;



% image plane
P.dx3 = P.lambda*P.z45/(P.Nx*P.dx4);
P.dy3 = P.lambda*P.z45/(P.Ny*P.dy4);
% 2D grid
[P.X3,P.Y3] = meshgrid(((1:P.Nx)-floor(P.Nx/2)-1)*P.dx3,((1:P.Ny)-floor(P.Ny/2)-1)*P.dy3);
X3_axis = ((1:P.Nx)-floor(P.Nx/2)-1)*P.dx3;
Y3_axis = ((1:P.Ny)-floor(P.Ny/2)-1)*P.dy3;

%  mirror and aperture plane
% pixel widths
P.dx2 = P.lambda*P.z34/(P.Nx*P.dx3);
P.dy2 = P.lambda*P.z34/(P.Ny*P.dy3);
% 2D grid
[P.X2,P.Y2] = meshgrid(((1:P.Nx)-floor(P.Nx/2)-1)*P.dx2,((1:P.Ny)-floor(P.Ny/2)-1)*P.dy2);
X2_axis = ((1:P.Nx)-floor(P.Nx/2)-1)*P.dx2;
Y2_axis = ((1:P.Ny)-floor(P.Ny/2)-1)*P.dy2;

% object plane
% pixel widths
P.dx1 = P.lambda*P.z12/(P.Nx*P.dx2);
P.dy1 = P.lambda*P.z12/(P.Ny*P.dy2);
% 2D grid
[P.X1,P.Y1] = meshgrid(((1:P.Nx)-floor(P.Nx/2)-1)*P.dx1,((1:P.Ny)-floor(P.Ny/2)-1)*P.dy1);
X1_axis = ((1:P.Nx)-floor(P.Nx/2)-1)*P.dx1;
Y1_axis = ((1:P.Ny)-floor(P.Ny/2)-1)*P.dy1;


P.P = 1/((1/P.f)-(1/(P.z12+P.z23))) ;
fprintf('Shifting the detector position\nchanges the resolution of the image plane\n,leading to a tighter or stretcher scan range\nin the 3rd section.\n')
%% Propagation & Image formation
% A point source model
P.D = 10e-5 ;
A = besselj(0,P.X1/P.D) .* besselj(0,P.Y1/P.D).* circ(P.X1,P.Y1,5*P.D) ;
% propagation of source radaition
Psi2 =  Afactor(P.z12,P.k,P.X2,P.Y2) .* fourier2D(A.*Bfactor(P.z12,P.k,P.X1,P.Y1));
% A circular aperture 
P.mask = circ(P.X2,P.Y2,P.O) ;
% Aperture's effect
Psi2 = Psi2 .* P.mask ;
% Propagation steps :
% 1 -  aperture to mirror
Psi3 =  fnear1D(P.Nx,P.dx2,P.k,P.z23,Psi2) ;
% 2 - mirror to image
Psi4 =  Afactor(P.z34,P.k,P.X3,P.Y3) .* fourier2D(Psi3.*Bfactor(-P.f,P.k,P.X2,P.Y2).*Bfactor(P.z34,P.k,P.X2,P.Y2));
% 3 - image to detector
Psi5 =  Afactor(P.z45,P.k,P.X4,P.Y4) .* fourier2D(Psi4.*Bfactor(P.z45,P.k,P.X3,P.Y3));

% Back propagation steps: 
% 1 -  detector to image
Psi6 =  -Afactor(-P.z45,P.k,P.X3,P.Y3) .* fourier2D(Psi5.*Bfactor(-P.z45,P.k,P.X4,P.Y4)) ;
% 2 -  image to mirror
Psi7 =  Afactor(-P.z34,P.k,P.X2,P.Y2)  .* fourier2D(Psi6.*Bfactor(-P.z34,P.k,P.X3,P.Y3)) ...
.* Bfactor(P.f,P.k,P.X2,P.Y2) ;
% 3 -  mirror to aperture 
Psi8 =  fnear1D(P.Nx,P.dx2,P.k,-P.z23,Psi7) .* P.mask ;
% Plot(Forward Propagation)
 figure;imagesc(A),axis equal tight,
 figure;imagesc(abs(Psi2).^2),axis equal tight,%colorbar,title('aperture')
 figure;imagesc(abs(Psi3).^2),axis equal tight,%colorbar,title('mirror')
 figure;imagesc(abs(Psi4).^2),axis equal tight,%colorbar,title('image')
 figure;imagesc(abs(Psi5).^2),axis equal tight,%colorbar,title('detector')
% Plot(Back Propagation)
 figure;imagesc(abs(Psi6).^2),axis equal tight,colorbar,%title('image')
 figure;imagesc(abs(Psi7).^2),axis equal tight,colorbar,%title('mirror')
 %figure;imagesc(abs(Psi8).^2),axis equal tight,colorbar,

 %% Initialization(an object suggestion)
% Emperically measured intensity
P.I_exp = abs(Psi5).^2 ;
% normalization factor
fnorm = sqrt(P.Nx*P.Ny);
% A random gaussian phase distribution
P.sigma = 5 ;
P.alpha = 1 ;
G = fspecial('gaussian',5*P.sigma,2*P.sigma);
phase_0 = imfilter(rand(P.Ny,P.Nx),G,'circular'); 
phase_0 = (phase_0-min(phase_0(:)))/(max(phase_0(:)-min(phase_0(:))));
% A primary source field suggestion
%P.s = 0.5e-4 ;
%Ki0 = exp(-(P.X1.^2+P.Y1.^2)/(2*P.s^2)) ;
% A bad guess
Ki0 = real(ifourier2D(ifourier2D(P.I_exp))) ;
Ki0 = Ki0 .* exp(1i*P.alpha*2*pi*phase_0) ; 
% Propgation to the aperture plane
%sf = 0.2 ;%source-position fluctuation
%d = P.z12 + sf ;
%d = 100 ;
%Ki0 = Afactor(d,P.k,P.X2,P.Y2) .* fourier2D(Ki0.*Bfactor(P.k,d,P.X1,P.Y1))/fnorm ;
% Aperture effect 
%Ki0 = Psi2 ;
%Ki0 = Ki0 .* P.mask ;
figure;imagesc(abs(Ki0).^2),axis equal tight
figure;imagesc(angle(Ki0)),axis equal tight

%% Reconstruction (from object plane)
for i = 1:1000
    disp(i)
Ki1 =  fnear1D(P.Nx,P.dx2,P.k,P.z23,Ki0) ;
Ki2 =  Afactor(P.z34,P.k,P.X3,P.Y3) .* fourier2D(Ki1);
Ki3 =  -1i*exp(1i*P.k*P.z45) * fourier2D(Ki2.*Bfactor(P.z45,P.k,P.X3,P.Y3));
Ki4 =  sqrt(P.I_exp).* exp(1i*angle(Ki3)) ;
Ki5 =  Afactor(-P.z45,P.k,P.X3,P.Y3) .* fourier2D(Ki4) ;
Ki6 =  -1i*exp(-1i*P.k*P.z34)  * fourier2D(Ki5.*Bfactor(-P.z34,P.k,P.X3,P.Y3)) ;
Ki7 =  fnear1D(P.Nx,P.dx2,P.k,-P.z23,Ki6) ;
Ki7 =  Ki7 .* P.mask  ;
Ki0 =  Ki7 ;

imagesc(angle(Ki0)),axis equal tight,drawnow

end
% figure;imagesc(abs(Ki1).^2),axis equal tight,colorbar
 figure;imagesc(abs(Ki2).^2),axis equal tight,colorbar
 figure;imagesc(abs(Ki3).^2),axis equal tight,colorbar
 
 
 %%

%  
% % Emperically measured intensity
% P.I_exp = abs(Psi5).^2 ;
% % normalization factor
% fnorm = sqrt(P.Nx*P.Ny);
% % A random gaussian phase distribution
% P.sigma = 5 ;
% P.alpha = 1 ;
% G = fspecial('gaussian',5*P.sigma,2*P.sigma);
% phase_0 = imfilter(rand(P.Ny,P.Nx),G,'circular'); 
% phase_0 = (phase_0-min(phase_0(:)))/(max(phase_0(:)-min(phase_0(:))));
% % A primary source field suggestion
% %P.s = 0.5e-4 ;
% %Ki0 = exp(-(P.X1.^2+P.Y1.^2)/(2*P.s^2)) ;
% % A bad guess
% Ki1 = abs(Psi5) ;
% Ki1 = Ki1 .* exp(1i*P.alpha*2*pi*phase_0); 
% %Ki1 = Psi5 .* conj(Afactor(P.z45,P.k,P.X4,P.Y4)) ;
% figure;imagesc(abs(Ki1).^2),axis equal tight
% figure;imagesc(angle(Ki1)),axis equal tight
%  %%
%  for i = 1:100
%     disp(i)
% Ki2 =   Afactor(-P.z45,P.k,P.X3,P.Y3) .* fourier2D(Ki1) ;
% Ki3 = -1i*exp(-1i*P.k*P.z34) .* fourier2D(Ki2.*Bfactor(-P.z34,P.k,P.X3,P.Y3)) ;
% Ki4 = fnear1D(P.Nx,P.dx2,P.k,-P.z23,Ki3) ;   
% Ki4 = Ki4 .* P.mask  ;   
% Ki5 =  fnear1D(P.Nx,P.dx2,P.k,P.z23,Ki4) ;
% Ki6 =  Afactor(P.z34,P.k,P.X3,P.Y3) .* fourier2D(Ki5);
% Ki7 =  -1i*exp(1i*P.k*P.z45) .* fourier2D(Ki6.*Bfactor(P.z45,P.k,P.X3,P.Y3));  
% Ki7 = sqrt(P.I_exp) .* exp(1i*angle(Ki7)) ;   
% Ki1 = Ki7 ;
% 
% imagesc(angle(Ki1)),axis equal tight,drawnow
% 
% end
% %figure;imagesc(abs(Ki1).^2),axis equal tight,colorbar
% %figure;imagesc(abs(Ki2).^2),axis equal tight,colorbar
% %figure;imagesc(abs(Ki3).^2),axis equal tight,colorbar
%  
%  
 
 
 