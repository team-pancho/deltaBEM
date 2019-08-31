function [beta0,beta1] = CalderonCalculusTest(u,gradu,gp,gm,varargin)

% [beta0,beta1] = CalderonCalculusTest(u,gradu,gp,gm)
% [beta0,beta1] = CalderonCalculusTest(u,gradu,gp,gm,fork)
% Input:
%       u : vectorized scalar function of x,y
%   gradu : vectorized vector-valued function of x,y (returns K x 2 matrix)
%   gp,gm : two companion meshes
%    fork : 0 or absent (averaged method) ~=0 (mixed method)
% Output:
%   beta0 : Dirichlet boundary data
%   beta1 : Neumann boundary data
%
% Last Update: October 31, 2013

% Tests on two meshes

beta0p = u(gp.midpt(:,1),gp.midpt(:,2));
beta1p = sum(gradu(gp.midpt(:,1),gp.midpt(:,2)).*gp.normal,2);  

beta0m = u(gm.midpt(:,1),gm.midpt(:,2));
beta1m = sum(gradu(gm.midpt(:,1),gm.midpt(:,2)).*gm.normal,2);  

% averaging/mixing and quadrature

if nargin==5
    fork=varargin{1};
else
    fork=0;
end
[Q,~,Pp,Pm] = CalderonCalculusMatrices(gp,fork);
beta0 = Pp*beta0p+Pm*beta0m;
beta1 = Q*(Pp*beta1p+Pm*beta1m);

return