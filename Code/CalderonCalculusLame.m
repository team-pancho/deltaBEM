function [V,K,J,W,C] = CalderonCalculusLame(g,gp,gm,mu,lambda)

% [V,K,J,W,C] = CalderonCalculusLame(g,gp,gm,mu,lambda)
% Input:
%     g : principal geometry
%    gp : companion geometry with epsilon = 1/6
%    gm : companion geometry with epsilon = -1/6
%    mu,lambda : Lame parameters
% Output:
%     V : single layer operator 
%     K : double layer operator K 
%     J : transpose of operator K 
%     W : hypersingular operator W
%     C : projection on local rigid motions
%
% Last Modified: November 4, 2013

[Vp,Kp,Jp,Wp,Cp] = CalderonCalculusLameHalf(g,gp,mu,lambda);
[Vm,Km,Jm,Wm,Cm] = CalderonCalculusLameHalf(g,gm,mu,lambda);

[Q,~,Pp,Pm] = CalderonCalculusMatrices(g,1);
O = zeros(size(Q));
Q  = [Q O; O Q]; 
Pp = [Pp O;O Pp];
Pm = [Pm O;O Pm];
V = Pp*Vp+Pm*Vm;
K = (Pp*Kp+Pm*Km)*Q;
J = Q*(Pp*Jp+Pm*Jm);
W = Pp*Wp+Pm*Wm;
C = Q*(Pp*Cp+Pm*Cm)*Q;

return

% Subfunction computing the two halves of the operators

function [V,K,J,W,P] = CalderonCalculusLameHalf(g,gp,mu,lambda)

% [V,K,J,W,P] = CalderonCalculusLameHalf(g,gp,mu,lambda)
% Input:
%     g : principal geometry
%    gp : companion geometry
%    mu,lambda : Lame parameters
% Output:
%     V : single layer operator
%     K : double layer operator
%     J : transpose of double layer operator
%     W : hypersingular operator 
%     P : projection on rigid motions
%
% Last Modified: November 4, 2013

A = mu/(lambda+2*mu);

% V and K = First row of Calderon Projector

RX  = bsxfun(@minus,gp.midpt(:,1),g.midpt(:,1)');
RY  = bsxfun(@minus,gp.midpt(:,2),g.midpt(:,2)');
RXNX = bsxfun(@times,RX,g.normal(:,1)');
RXNY = bsxfun(@times,RX,g.normal(:,2)');
RYNX = bsxfun(@times,RY,g.normal(:,1)');
RYNY = bsxfun(@times,RY,g.normal(:,2)');
RN   = RXNX+RYNY;
RT   = -RXNY+RYNX;
R = sqrt(RX.^2+RY.^2);      
O = zeros(size(R));

V = (1+A)/(4*pi*mu)*[-log(R) O; O -log(R)]...
    +(1-A)/(4*pi*mu)*repmat(1./R.^2,[2 2])...
                     .*[RX.*RX  RX.*RY; RY.*RX RY.*RY];
K = A/(2*pi)*repmat(1./R.^2,[2 2]).*([RN O; O RN]+[O RT; -RT O])...
    +(1-A)/pi*repmat(RN./R.^4,[2 2]).*[RX.*RX  RX.*RY; RY.*RX  RY.*RY];

% J = Transposed double layer operator

RX = -RX; RY = -RY;
NXRX = bsxfun(@times,gp.normal(:,1),RX);
NYRX = bsxfun(@times,gp.normal(:,2),RX);
NXRY = bsxfun(@times,gp.normal(:,1),RY);
NYRY = bsxfun(@times,gp.normal(:,2),RY);
NR   = NXRX+NYRY;
TR   = -NYRX+NXRY;

J = A/(2*pi)*repmat(1./R.^2,[2 2]).*([NR O; O NR]+[O -TR; TR O])...
    +(1-A)/pi*repmat(NR./R.^4,[2 2]).*[RX.*RX  RX.*RY; RY.*RX  RY.*RY];

% W = Hypersingular operator

RX = bsxfun(@minus,gp.brkpt(:,1),g.brkpt(:,1)');
RY = bsxfun(@minus,gp.brkpt(:,2),g.brkpt(:,2)');
R  = sqrt(RX.^2+RY.^2);                          % |b_i^ep-b_j|
W  = A*(lambda+mu)/pi*([-log(R) O;O -log(R)]...
         +repmat(1./R.^2,[2 2]).*[RX.*RX  RX.*RY; RY.*RX  RY.*RY]);
next = [g.next length(g.next)+g.next];   % next index in two components
nextp = [gp.next length(gp.next)+gp.next];  

W = W(nextp,next)-W(nextp,:)-W(:,next)+W;
            
% Projection onto rigid motions

l=sum(g.normal.^2,2); lp=sum(gp.normal.^2,2);
x=l.*g.midpt(:,1); xp=lp.*gp.midpt(:,1);
y=l.*g.midpt(:,2); yp=lp.*gp.midpt(:,2);
N = size(g.midpt,1);
g.comp=[g.comp N+1];
C = sparse(N,N);
XX = C; XY = C; YX = C; YY= C;
for c=1:length(g.comp)-1
    list=g.comp(c):g.comp(c+1)-1;
    C(list,list)=lp(list)*l(list)';
    XX(list,list)=xp(list)*x(list)';
    YY(list,list)=yp(list)*y(list)';
    XY(list,list)=xp(list)*y(list)';
    YX(list,list)=yp(list)*x(list)';
end
O = sparse(N,N);
P = [C O; O C]+[YY -YX; -XY XX];  

return


