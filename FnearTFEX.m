function f1 = FnearTFEX(dx,nx,lambda,z,f0)

 zcrit = nx * dx^2/lambda ;
 if abs(z) < zcrit
     
     k = 2*pi/lambda ;
     du = 2*pi / (nx * dx) ;
     [ux,uy] = meshgrid(du *(-nx/2:nx/2-1),du * (-nx/2:nx/2-1)) ;
     H = exp(-1i*(ux.^2+uy.^2)*z/(2*k));
     f1 = ifourier2D(fourier2D(f0) .* H);
     f1 = exp(1i*k*z) * f1 ;
    
 else
     fprintf('zcirt = %d \n',zcrit)
 end
 
