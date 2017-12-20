

function [sharpness,z_image,E_c,Scan] = caustic(dx,M,lambda,E) 
 
  Range = 1.9 * dx^2 .*  M/lambda ;
  NN = 250;
  dz = Range/NN ;

  z_range = ((1:NN)-floor(NN/2)-1)*dz;
  PROBE = zeros(M,M,NN);
  sharpness = zeros(NN,1);

   for jj = 1:NN
      disp(jj) ;
      PROBE(:,:,jj) = FnearTFEX(dx,M,lambda,z_range(jj),E);
      P4 = abs((PROBE(:,:,jj))).^4;
      sharpness(jj) = sum(P4(:));
   end
   
   Scan = reshape(PROBE(M/2,:,:),M,NN);
   Scan = abs(Scan).^2 ;
   Scan = Scan/max(Scan(:)) ;
   z_image = z_range(sharpness == max(sharpness));
   E_c = FnearTFEX(dx,M,lambda,z_image,E);

end