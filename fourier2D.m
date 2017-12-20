% ftt = fourier2D(f) is a 2D Fourier transfrom.
% INPUT:
% f : an arbitrary function
% OUTPUT:
% ftt : 2D matrix specifying the Fourier transform of an arbitrary function.

% M. Mehrjoo, K. Giewekemeyer, European XFEL (2015)

function ftt = fourier2D(f)

ftt = fftshift(gather(fft2(ifftshift(f))));