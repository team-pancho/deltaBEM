function [Q,M,Pp,Pm] = CalderonCalculusMatrices(g,varargin)

% [Q,M,Pp,Pm] = CalderonCalculusMatrices(g)
% [Q,M,Pp,Pm] = CalderonCalculusMatrices(g,fork)
% Input : 
%       g : geometry
%    fork : 0 or absent (averaged method) ~=0 (mixed method)
% Output:
%      Q  : lookaround quadrature matrix
%      M  : mass matrix
%   Pp,Pm : mixing matrices (=1/2 if fork=0)
%
% Last Modified: October 31, 2013

fork=0;
if nargin==2
    fork=varargin{1};
end
if fork 
    alpha=5/6;
else
    alpha=1;
end

N = size(g.midpt,1);

% quadrature matrix

Q = sparse(1:N,1:N,22/24)...
    +sparse(g.next,1:N,1/24)+sparse(1:N,g.next,1/24);

% mass matrix

M = sparse(1:N,1:N,(4+12*alpha)/18)...
    +sparse(g.next,1:N,(7-6*alpha)/18)+sparse(1:N,g.next,(7-6*alpha)/18);

% mixing matrices or average values

if fork
    Pp=sparse(1:N,1:N,alpha/2)+sparse(g.next,1:N,(1-alpha)/2);
    Pm=Pp';
else
    Pp=1/2;
    Pm=1/2;
end

return

% singular quadrature

S = sparse(1:N,1:N,349/432)...
    +sparse(1:N,g.next,5/1296)...
    +sparse(g.next,1:N,89/432)...
    +sparse(g.next(g.next),1:N,-23/1296);
