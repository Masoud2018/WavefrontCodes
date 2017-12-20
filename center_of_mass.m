function [cx,cy] = center_of_mass(I)

% Author K. Giewekemeyer
%
% Last modified 10.12.2010

[x,y] = meshgrid(1:size(I,2),1:size(I,1));
cx = sum(sum(x.*I))/sum(sum(I));
cy = sum(sum(y.*I))/sum(sum(I));