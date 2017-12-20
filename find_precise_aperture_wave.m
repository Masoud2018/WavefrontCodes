

function f = find_precise_aperture_wave(dx1,M,lambda,E,ap)

%%%%%% Aperture's precise position %%%%%%

  Range = 1.9 * dx1^2 .* M / lambda ;
  NN = 100;
  dz = Range/NN ;

  z_range = ((1:NN)-floor(NN/2)-1)*dz;
  PROBE = zeros(M,M,NN);

   for jj = 1:NN
      disp(jj) ;
      PROBE(:,:,jj) = FnearTFEX(dx1,M,lambda,z_range(jj),E) ;
   end

fff = PROBE(:,:,(NN/2)-1) ;
ff = abs(fff).^2 .* ap ;
ff = smooth2D(ff,2,1) .* ap ;
f = sqrt(ff) .* exp(1i*angle(fff)) ;