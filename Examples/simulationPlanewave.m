% A script to make a simulation of transient scattering
% Last modified: December 22, 2015

% discretization parameters

N = 250;
M = 300;
T = 5;
Tlag = 1;
fork = 0;

% Trapezoid rule method transfer function & time step
p =@(z) 2*(1-z)/(1+z);
kappa = T/M;

% signal
HH =  @(x) x.^5.*(1-5*(x-1)+15*(x-1).^2-35*(x-1).^3+70*(x-1).^4-...
           126*(x-1).^5).*(x>0).*(x<1)+(x>=1);
HHp = @(x) (5*x.^4.*...
           (1-5*(x-1)+15*(x-1).^2-35*(x-1).^3+70*(x-1).^4-126*(x-1).^5)...
           +x.^5.*(-5+30*(x-1)-105*(x-1).^2+280*(x-1).^3-630*(x-1).^4))...
           .*(x>0).*(x<1);
alpha=1/8; l=3/8;   % [0,alpha] - [alpha,l-alpha] [l-alpha l]
window  = @(x) HH(x/alpha).*HH((l-x)./alpha);
windowp = @(x) 1/alpha*(HHp(x/alpha).*HH((l-x)./alpha)...
                        -HH(x/alpha).*HHp((l-x)./alpha));
signal = @(t) sin(16*t).*window(t); 
signalp= @(t) sin(16*t).*windowp(t)+16*cos(16*t).*window(t);

% direction of wave propogation
dir = [1 0];

% Discrete geometry
g  = affine(kite(N,0),0.5*eye(2),[0;0]);
gp = affine(kite(N,1/6),0.5*eye(2),[0;0]);
gm = affine(kite(N,-1/6),0.5*eye(2),[0;0]);

% Generate the signal and boundary data
[~,~,Pp,Pm] = CalderonCalculusMatrices(g,fork);
[uincp,dnuincp] = planewave(signal,signalp,dir,Tlag,gp.midpt,gp.normal,T,M,1);
[uincm,dnuincm] = planewave(signal,signalp,dir,Tlag,gm.midpt,gm.normal,T,M,1);
beta0 = (Pp*uincp+Pm*uincm);

% Solve and postprocess
[V,~,~,~] = CalderonCalculusHelmholtz(g,gp,gm,fork);
eta = -CQequation(V,beta0,kappa,p);

% plotting parameters
box = [-2 2 2 -2];
dabs = 0.1;
drel = 0.1;
hgrid = 0.03;   % mesh size
fill = 0;   %plot outside only
[X,Y,Tri,~,~]=triangulateGeometry(box,g,dabs,drel,hgrid,fill);

% domain decomp trick
Block=200;
Nvertices = size(X,1);
Npieces   = ceil(Nvertices/Block)
pieces = @(n) 1+Block*(n-1):min(Block*n,Nvertices);
uscat = zeros(size(X,1),M+1);

for n=1:Npieces
    disp(n)
    XX=X(pieces(n),:);
    YY=Y(pieces(n),:);
    [S,~]=HelmholtzPotentials(g,[XX YY]);
    uscat(pieces(n),:) = CQforward(S,eta,kappa,p);
end

uinc = planewave(signal,signalp,dir,Tlag,[X Y],0*[X Y],T,M,0);
utotal = uinc+uscat;

% pin the corners
TRC = find(X==max(X) & Y==max(Y));
BLC = find(X==min(X) & Y==min(Y));
utotal(TRC,:) = max(utotal(:));
utotal(BLC,:) = min(utotal(:));

for t=1:M+1
    trisurf(Tri(:,1:3),X,Y,utotal(:,t)), axis off, axis equal, view(2),
    shading interp, shg
    pause(0.05)
end
