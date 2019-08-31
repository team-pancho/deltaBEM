% A script to make a simulation of transient scattering
% Last modified: December 22, 2015

% discretization parameters

N = 150;
M = 350;
T = 5;
fork = 0;

% Trapezoid rule method transfer function
p =@(z) 2*(1-z)/(1+z);
kappa = T/M;

% signal
signal =@(t) sin(16*t).^5.*(t>=0);

% Discrete geometry
g1 = ellipse(N,0,[1 2]/2,[1 1]);
g1p = ellipse(N,1/6,[1 2]/2,[1 1]);
g1m = ellipse(N,-1/6,[1 2]/2,[1 1]);

g2 = ellipse(N,0,[2 1]/2,-[1 1]);
g2p = ellipse(N,1/6,[2 1]/2,-[1 1]);
g2m = ellipse(N,-1/6,[2 1]/2,-[1 1]);

g  = joinGeometry(g1,g2);
gp = joinGeometry(g1p,g2p);
gm = joinGeometry(g1m,g2m);

% Generate the signal and boundary data
[Q,~,Pp,Pm] = CalderonCalculusMatrices(g,fork);
src = [-1 0.65];
[~,~,~,W] = CalderonCalculusHelmholtz(g,gp,gm,fork);
[~,dnuincp]=cylindricalwave(signal,1,src,gp.midpt,gp.normal,T,M,1);
[~,dnuincm]=cylindricalwave(signal,1,src,gm.midpt,gm.normal,T,M,1);
beta1 = (Pp*dnuincp+Pm*dnuincm);

% Solve and postprocess
psi = CQequation(W,beta1,kappa,p);

% plotting parameters
box = [-2.5 2.5 2.5 -2.5];
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
    [~,D]=HelmholtzPotentials(g,[XX YY]);
    uscat(pieces(n),:) = CQforward(D,Q*psi,kappa,p);
end

uinc = cylindricalwave(signal,1,src,[X Y],0*[X Y],T,M,0);
utotal = uinc+uscat;

% pin the corners
TRC = find(X==max(X) & Y==max(Y));
BLC = find(X==min(X) & Y==min(Y));
utotal(TRC,:) = max(utotal(:));
utotal(BLC,:) = min(utotal(:));

for t=1:M+1
    trisurf(Tri(:,1:3),X,Y,utotal(:,t)), axis off, axis equal, view(2),
    shading interp, colormap(bone), shg
    pause(0.05)
end

%%

name='cylindricalMovie/frame';

Nindices=M+1;
indices=cell(1,Nindices);

for i=1:Nindices
    indices{i}=num2str(i,'%03d');
end


for i=1:length(indices)
    longname{i}=[name,indices{i},'.png'];
end
for t=1:M
    trisurf(Tri(:,1:3),X,Y,utotal(:,t+1)), shading interp, axis off, ...
        view(2), colormap(bone), axis equal, shg;
    pause(0.2)
    saveas(gcf,longname{t})
end
