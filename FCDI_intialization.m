
function f = FCDI_intialization(data,sigma,alpha,N)

I_exp = data ;
G = fspecial('gaussian',5*sigma,2*sigma);
phase_0 = imfilter(rand(N,N),G,'circular'); 
phase_0 = (phase_0-min(phase_0(:)))/(max(phase_0(:)-min(phase_0(:))));
K0 = sqrt(I_exp) ; 
f = K0 .* exp(1i*alpha*phase_0) ;
