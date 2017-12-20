% B = AFACTOR(Z,K,X,Y) produces a numerically relevant factor for the
% iterative reconstruction algorithm for wavefield reconstrution according
% to Quiney et al. Nat. Phys. 2, 101 (2006).
%
% INPUT:
% Z : Propagation distance (in m)
% K : Wavenumber, 2*pi/lambda (in inverse meters)
% X,Y: Coordinate mesh in plane perpendicular to optical axis (in meters)
% OUTPUT:
% A : 2D matrix specifying the factor B.
%
% M. Mehrjoo, K. Giewekemeyer, European XFEL (2015)


function B = Bfactor(z,k,x,y)

B = exp(1i*0.5*k*(x.^2+y.^2)/z) ;