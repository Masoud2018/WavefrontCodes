
function f = Airy(x,a,lambda,z,d) 

Z = (pi/lambda) * a * (x-d)/z ;
f = (pi/lambda)*(a^2/z)* (besselj(1,Z)./Z) ;

end

