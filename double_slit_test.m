

% Young's double pinhole experiment simulation

%% Parameters are from Andrej's thesis

lambda = 1e-9; % wavelenght
w = 1e-6; % slit width
z = 1 ; % pinholes to detector
gamma = 0.49 ; % degree of spatial coherence
d = 5e-6; % slit separation
I0 = 1 ; % intensity at each pinhole (This value is subject to change.)
dx = 1e-3 ; % sampling 
n = 1e4 ; % number of samplinf point
x = dx * linspace(-2,2,n) ; % coordinate
R = 1/dx ; % rescale factor
%% Intensity simulation 

theta = ((pi/lambda) * d^2/z) - ((2*pi/lambda) * (d * x/z)) ; 
Partial_Coherence_Factor = 1 + gamma * cos(theta) ;
I_1 = I0 * Airy(x,w,lambda,z,d) ; 
I_D_0 = 2 * (pi/lambda) * (w^2/z) * I_1.^2 ;
I_D = I_D_0 .* Partial_Coherence_Factor ;
I_D_Normalized  = I_D/max(I_D) ;

figure,plot(1e3*x,I_D,'r')

%% Fitting 

% I_D_0 is a known value since the geometry and wavelength of the
% experiment are known.

for i = 1 : 100
    
    alpha = 0.01 * i ; % covers the reange (0,1) ! 
    
    Partial_Coherence_Factor_fit = 1 + alpha * cos(theta) ;
    I_D_fit = I_D_0 .* Partial_Coherence_Factor_fit ;
    I_D_fit_Normalized(i,:)  = I_D_fit/max(I_D_fit) ;
    
    % Error of fit 
    chi = I_D_Normalized - I_D_fit_Normalized(i,:) ;
    Error(i) = abs(sum(chi)/100) ;
    
end
 

for ii = 1 : 100
    
    gamma_fit(ii) = 0.01 * ii ; 
    
    if Error(ii) <= 1e-2 % indicates the best fit
        
        fprintf('gamma_fit:%2.2f\n',gamma_fit(ii))
        figure, plot(R*x,I_D_Normalized,'r'),hold on,plot(R*x(1:80:end),I_D_fit_Normalized(ii,1:80:end),'bo'), hold off 
        title(['Young"s double pinhole experiment, \gamma_{12}(0) = ',num2str(gamma)]), xlabel('x(mm)'), ylabel('I(x)(a.u)')
        legend('Simulation','Fit')
        dim = [.71 .2 .95 .6];
        str = ['\gamma_{12}^{fit}(0) = ',num2str(gamma_fit(ii))] ;
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        
    end
end



