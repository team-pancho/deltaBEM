function [V,K,J,W] = CalderonCalculusHelmholtz(g,gp,gm,varargin)

% [V,K,J,W] = CalderonCalculusHelmholtz(g,gp,gm)
% [V,K,J,W] = CalderonCalculusHelmholtz(g,gp,gm,fork)
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
%
% Last Modified: October 31, 2013

[Vp,Kp,Jp,Wpp,Vnp] = CalderonCalculusHelmholtzHalf(g,gp);
[Vm,Km,Jm,Wpm,Vnm] = CalderonCalculusHelmholtzHalf(g,gm);

fork=0;
if nargin==4
    fork=varargin{1};
end
[Q,~,Pp,Pm] = CalderonCalculusMatrices(g,fork);

V =@(s) Pp*Vp(s) + Pm*Vm(s);
K =@(s) (Pp*Kp(s) + Pm*Km(s))*Q;
J =@(s) Q*(Pp*Jp(s) + Pm*Jm(s));
W =@(s) Pp*Wpp(s) + Pm*Wpm(s) + Q*(Pp*Vnp(s) + Pm*Vnm(s))*Q;

return

% Subfunction computing the two halves of the Calculus

function [V,K,J,Wp,Vn] = CalderonCalculusHelmholtzHalf(g,gp)

% [V,K,J,Wp,Vn] = CalderonCalculusHelmholtzHalf(g,gp)
% Input:
%     g : principal geometry
%    gp : companion geometry
% Output:
%     V : single layer operator
%     K : double layer operator
%     J : transpose of double layer operator
%    Wp : principal part of hypersingular operator W
%    Vn : regular part of W
%
% Last Modified: September 5, 2013

% V(s)  = Single layer operator
% Vn(s) = Regular part of W(s)

DX = bsxfun(@minus,gp.midpt(:,1),g.midpt(:,1)');
DY = bsxfun(@minus,gp.midpt(:,2),g.midpt(:,2)');
D = sqrt(DX.^2+DY.^2);                    % |m_i^ep-m_j|
H  = gp.normal(:,1)*g.normal(:,1)'...
     +gp.normal(:,2)*g.normal(:,2)'; % n_i . n_j
 
V = @(s) 1i/4.*besselh(0,1,1i*s*D);
Vn = @(s) s^2.*H.*1i/4.*besselh(0,1,1i*s*D);

% Principal part of W(s)

DX = bsxfun(@minus,gp.brkpt(:,1),g.brkpt(:,1)');
DY = bsxfun(@minus,gp.brkpt(:,2),g.brkpt(:,2)');
D1 = sqrt(DX.^2+DY.^2);                          % |b_i^ep-b_j|
D2 = D1(gp.next,:);                             % |b_i+1^ep-b_j|   
D3 = D1(:,g.next);                               % |b_i^ep-b_j+1|
D4 = D1(gp.next,g.next);                        % |b_i+1^ep-b_j+1|
Wp = @(s) 1i/4*besselh(0,1,1i*s*D1)-1i/4*besselh(0,1,1i*s*D2)...
          -1i/4*besselh(0,1,1i*s*D3)+1i/4*besselh(0,1,1i*s*D4);

% K(s) = Double layer operator

DX = bsxfun(@minus,gp.midpt(:,1),g.midpt(:,1)');
DY = bsxfun(@minus,gp.midpt(:,2),g.midpt(:,2)');
D  = sqrt(DX.^2+DY.^2);           
N = bsxfun(@times,DX,g.normal(:,1)')...
    +bsxfun(@times,DY,g.normal(:,2)');
N = N./D;
K = @(s) -s/4.*besselh(1,1,1i*s*D).*N;

% J(s) = Transposed double layer operator
    
N = bsxfun(@times,gp.normal(:,1),DX)...
     +bsxfun(@times,gp.normal(:,2),DY);
N = N./D;
J = @(s) s/4.*besselh(1,1,1i*s*D).*N;

return


