function f = FCDI_intialization_inverse(data,sigma,alpha,N)

G = fspecial('gaussian',5*sigma,2*sigma);
phase_0 = imfilter(rand(N,N),G,'circular'); 
phase_0 = (phase_0-min(phase_0(:)))/(max(phase_0(:)-min(phase_0(:))));
f = abs(data) .* exp(1i*alpha*phase_0) ;