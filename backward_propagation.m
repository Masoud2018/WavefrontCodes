
function H = backward_propagation(F,I_exp,lambda,k,z34,z23,X2,Y2,X3,Y3,f)


g = (-1i/(lambda*z23)) * exp(1i*k*z23) * Bfactor(z23,k,X3,Y3) ...
    .* fourier2D(F .* Bfactor(-(z23*f)/(z23-f),k,X2,Y2)) ;
g = (-1i/(lambda*z34)) * exp(1i*k*z34)...
    .* fourier2D(g .* Bfactor(z34,k,X3,Y3)) ;
H =  sqrt(I_exp).* exp(1i*(angle(g))) ;