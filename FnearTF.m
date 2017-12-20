function f1 = FnearTF(L,dx,nx,lambda,z,f0)

 zcrit = L * dx/lambda ;
 if abs(z) < zcrit
     
     k = 2*pi/lambda ;
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