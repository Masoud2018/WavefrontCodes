
% c = circ(x,y,d) produces a circle with diameter D
%
% INPUT :
% X(Y) : Coordinate mesh 
% D : Circle diameter
% OUTPUT :
% A circle!



function c = circ(x,y,d)
r=sqrt(x.^2+y.^2);
c=double(r<d/2);
c(r==d/2)=1;

