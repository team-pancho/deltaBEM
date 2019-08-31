function [V,K,J,W,C] = CalderonCalculusLaplace(g,gp,gm,varargin)

% [V,K,J,W,C] = CalderonCalculusLaplace(g,gp,gm)
% [V,K,J,W,C] = CalderonCalculusLaplace(g,gp,gm,fork)
% Input:
%     g : principal geometry
%    gp : companion geometry with epsilon = 1/6
%    gm : companion geometry with epsilon = -1/6
%    fork : 0 or absent (averaged method) ~=0 (mixed method)
% Output:
%     V : single layer operator V
%     K : double layer operator K
%     J : transpose of operator K
%     W : hypersingular operator W
%     C : rank Ncomp perturbation
%
% Last Modified: October 31, 2013

[Vp,Kp,Jp,Wp,Cp] = CalderonCalculusLaplaceHalf(g,gp);
[Vm,Km,Jm,Wm,Cm] = CalderonCalculusLaplaceHalf(g,gm);

fork=0;
if nargin==2
    fork=varargin{1};
end
[Q,~,Pp,Pm] = CalderonCalculusMatrices(g,fork);
V = Pp*Vp+Pm*Vm;
K = (Pp*Kp+Pm*Km)*Q;
J = Q*(Pp*Jp+Pm*Jm);
W = Pp*Wp+Pm*Wm;
C = Q*(Pp*Cp+Pm*Cm)*Q;

return

% Subfunction computing the two halves of the operators

function [V,K,J,W,C] = CalderonCalculusLaplaceHalf(g,gp)

% [V,K,J,W,C] = CalderonCalculusLaplaceHalf(g,gp)
% Input:
%     g : principal geometry
%    gp : companion geometry
% Output:
%     V : single layer operator
%     K : double layer operator
%     J : transpose of double layer operator
%     W : hypersingular operator 
%     C : rank Ncomp perturbation
%
% Last Modified: September 4, 2013

% V = Single layer operator

DX = bsxfun(@minus,gp.midpt(:,1),g.midpt(:,1)');
DY = bsxfun(@minus,gp.midpt(:,2),g.midpt(:,2)');
D = sqrt(DX.^2+DY.^2);                    % |m_i^ep-m_j|
V = -1/(2*pi)*log(D);

% K = Double layer operator
      
N = bsxfun(@times,DX,g.normal(:,1)')...
    +bsxfun(@times,DY,g.normal(:,2)');
K = 1/(2*pi)*N./D.^2;

% J = Transposed double layer operator
    
N = bsxfun(@times,gp.normal(:,1),DX)...
     +bsxfun(@times,gp.normal(:,2),DY);
J = -1/(2*pi)*N./D.^2;

% W = Hypersingular operator

DX = bsxfun(@minus,gp.brkpt(:,1),g.brkpt(:,1)');
DY = bsxfun(@minus,gp.brkpt(:,2),g.brkpt(:,2)');
D  = sqrt(DX.^2+DY.^2);                          % |b_i^ep-b_j|
W  = -1/(2*pi)*(log(D(gp.next,g.next))+log(D)...
                -log(D(gp.next,:))-log(D(:,g.next)));
            
% C = rank Ncomp perturbation

lengths=sum(g.normal.^2,2);
lengthsp=sum(gp.normal.^2,2);
N = size(g.midpt,1);
g.comp=[g.comp N+1];
C =sparse(N,N);
for c=1:length(g.comp)-1
    list=g.comp(c):g.comp(c+1)-1;
    C(list,list)=lengthsp(list)*lengths(list)';
end

return


