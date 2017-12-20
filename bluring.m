

function [f,Gamma] = bluring(NN,sigma,I)


[X,Y] = meshgrid(((1:NN)-floor(NN/2)-1),((1:NN)-floor(NN/2)-1));
Gamma = exp(-(X.^2+Y.^2)/(2*sigma.^2)) ;
f = conv2(I,Gamma,'same') ;
f = norm_function(f) ;
