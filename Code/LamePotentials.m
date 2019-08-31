function [SL,DL]=LamePotentials(g,z,mu,lambda)

% [SL,DL]=LamePotentials(g,z,mu,lambda)
% 
% Input:
%     g : geometry
%     z : K x 2 matrix with points where potentials are evaluated
%     mu, lambda : Lame parameters
% Output:
%   SL  : matrix for SL
%   DL  : matrix for DL
% Last modified: November 5, 2013

A = mu/(lambda+2*mu);

RX  = bsxfun(@minus,z(:,1),g.midpt(:,1)');
RY  = bsxfun(@minus,z(:,2),g.midpt(:,2)');
RXNX = bsxfun(@times,RX,g.normal(:,1)');
RXNY = bsxfun(@times,RX,g.normal(:,2)');
RYNX = bsxfun(@times,RY,g.normal(:,1)');
RYNY = bsxfun(@times,RY,g.normal(:,2)');
RN   = RXNX+RYNY;
RT   = -RXNY+RYNX;
R = sqrt(RX.^2+RY.^2);  
O = zeros(size(R));

SL = (1+A)/(4*pi*mu)*[-log(R) O; O -log(R)]...
     +(1-A)/(4*pi*mu)*repmat(1./R.^2,[2 2])...
                     .*[RX.*RX  RX.*RY; RY.*RX RY.*RY];                  
DL = A/(2*pi)*repmat(1./R.^2,[2 2]).*([RN O; O RN]+[O RT; -RT O])...
     +(1-A)/pi*repmat(RN./R.^4,[2 2])...
              .*[RX.*RX RX.*RY; RY.*RX RY.*RY];

return