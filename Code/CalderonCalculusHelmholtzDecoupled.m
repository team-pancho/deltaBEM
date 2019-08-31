function [V,K,J,W]=CalderonCalculusHelmholtzDecoupled(g,gp,gm,c,alpha,varargin)

% [V,K,J,W]=CalderonCalculusHelmholtzDecoupled(g,gp,gm,c,alpha)
% [V,K,J,W]=CalderonCalculusHelmholtzDecoupled(g,gp,gm,c,alpha,fork)
%
% Input:
%    g,gp,gm   : geometry structures for main and companion meshes
%                for a geometry with K components
%    c, alpha  : K x 1 vectors
%    fork      : 0 or absent (averaged method) ~=0 (mixed method)
% Output:
%    V,K,J,W   : block diagonal operators (Costabel-Stephan formulation)
%                chi^{-1} V(s/c), K(s/c), J(s/c), chi W(s/c)
% Last modified: October 31, 2013

G=unpackGeometry(g);
Gp=unpackGeometry(gp);
Gm=unpackGeometry(gm);
nComp=length(G);

fork=0;
if nargin==6
    fork=varargin{1};
end
for comp=1:nComp
    [VV,KK,JJ,WW]=CalderonCalculusHelmholtz(G{comp},Gp{comp},Gm{comp},fork);
     V{comp}=VV;
     K{comp}=KK;
     J{comp}=JJ;
     W{comp}=WW;
end
V=@(s) attachBlockDiag(V,s,c,alpha.^(-1));
W=@(s) attachBlockDiag(W,s,c,alpha);
J=@(s) attachBlockDiag(J,s,c,ones(size(alpha)));
K=@(s) attachBlockDiag(K,s,c,ones(size(alpha)));
return

% Subfunction to merge functions in a diagonal block function

function F=attachBlockDiag(f,s,speed,factor)
a=cell(1,length(f));
for i=1:length(f);
    a{i}=factor(i)*f{i}(s/speed(i));
end
F=blkdiag(a{:});
return
