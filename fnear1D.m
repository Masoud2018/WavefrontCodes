function f1 = fnear1D(nx,dx,k,z,f0)

 zcrit= nx * dx^2 /(2 * pi / k) ;
 %zcrit = L^2/(nx*(2 * pi / k)) ;
 if abs(z) < zcrit
                                    
     du = 2*pi / (nx * dx) ;
     [ux,uy] = meshgrid(du *(-nx/2:nx/2-1),du * (-nx/2:nx/2-1)) ;
     H=exp(-1i*(ux.^2+uy.^2)*z/(2*k));
     f1 = ifourier2D(fourier2D(f0) .* H);
     f1 = exp(1i*k*z) * f1 ;
    
    
 else
     fprintf('zcirt = %d \n',zcrit)
 end
 
 function ftt = fourier2D(f)

ftt = fftshift(fft2(ifftshift(f)));

function ftt = ifourier2D(f)

ftt = fftshift(ifft2(ifftshift(f)));