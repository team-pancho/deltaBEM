function [SL,DL]=HelmholtzPotentials(g,z)

% [SL,DL]=HelmholtzPotentials(g,z)
% 
% Input:
%     g : geometry
%     z : K x 2 matrix with points where potentials are evaluated
% Output:
%   SL(s)     : transfer function for SL
%   DL(s)     : transfer function for DL
% Last modified: May 28, 2014.

% SL(s) : single layer potential

RX = bsxfun(@minus,z(:,1),g.midpt(:,1)');   
RY = bsxfun(@minus,z(:,2),g.midpt(:,2)');
R = sqrt(RX.^2+RY.^2);  
SL = @(s) 1i/4*besselh(0,1,1i*s*R);

% DL(s) : double layer potential

RN = bsxfun(@times,RX,g.normal(:,1)')...
    +bsxfun(@times,RY,g.normal(:,2)');
RN = RN./R;
DL = @(s) -s/4*besselh(1,1,1i*s*R).*RN;

return