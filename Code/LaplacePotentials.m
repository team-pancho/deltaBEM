function [SL,DL]=LaplacePotentials(g,z)

% [SL,DL]=LaplacePotentials(g,z)
% 
% Input:
%     g : geometry
%     z : K x 2 matrix with points where potentials are evaluated
% Output:
%   SL  : matrix for SL
%   DL  : matrix for DL
% Last modified: September 4, 2013

% SL : single layer potential

DX = bsxfun(@minus,z(:,1),g.midpt(:,1)');   
DY = bsxfun(@minus,z(:,2),g.midpt(:,2)');
D = sqrt(DX.^2+DY.^2);  
SL = -1/(2*pi)*log(D);

% DL : double layer potential

N = bsxfun(@times,DX,g.normal(:,1)')...
    +bsxfun(@times,DY,g.normal(:,2)');
DL = 1/(2*pi)*N./D.^2;

return