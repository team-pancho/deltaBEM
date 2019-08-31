% Script for a scattering problem by two obstacles, with 
% Dirichlet and Neumann boundary conditions

% Last modified: January 12, 2015

if ~exist('external')
    k = input('Input frequency: ');
    N = input('Number of points: ');
end

% FREQUENCY

s=-1i*k;

% EXACT SOLUTION

sc  = [0.1 0.2];    % source point, interior to one obstacle
r   = @(x,y) sqrt((x-sc(1)).^2+(y-sc(2)).^2);
u   = @(x,y) besselh(0,1,1i*s*r(x,y));    % Exterior radiating solution
ux  = @(x,y) -(x-sc(1)).*1i*s.*besselh(1,1,1i*s*r(x,y))./r(x,y);
uy  = @(x,y) -(y-sc(2)).*1i*s.*besselh(1,1,1i*s*r(x,y))./r(x,y);
gradu = @(x,y) [ux(x,y) uy(x,y)];

% INCIDENT WAVE

uinc = @(x,y) -u(x,y);
graduinc = @(x,y) -gradu(x,y);

% DISCRETE GEOMETRY
% Sound-soft ellipse, tagged with index 1
% Sound-hard kite, tagged with index 2

[g1,g1p,g1m]=sample(@ellipse,N,[1.2 1],[3,0]);

[g2,g2p,g2m]=sample(@kite,N);
[g2,g2p,g2m]=threetimes(@affine,g2,g2p,g2m,0.5*eye(2),[0;0]);

[g,gp,gm]=merge({g1,g2},{g1p,g2p},{g1m,g2m});

% Restriction Matrices

[g1indices,R1] = selectComponents(g,1);
[g2indices,R2] = selectComponents(g,2);

% OBSERVATION POINTS

xobs=[1.5 0;-1 0;1 -1];

% OPERATORS, POTENTIALS, AND RIGHT-HAND SIDES

[V1,~,~,~] = CalderonCalculusHelmholtz(g1,g1p,g1m);
[~,~,~,W2] = CalderonCalculusHelmholtz(g2,g2p,g2m);
if isfield(g1,'parity')   % In case g1 is an open arc
    V1 = @(s) V1(s) + g1.parity;
end
if isfield(g2,'parity')   % In case g2 is an open arc
    W2 = @(s) W2(s) + abs(g2.parity);
end
Q2 = CalderonCalculusMatrices(g2);
[~,K,J,~] = CalderonCalculusHelmholtz(g,gp,gm);
[beta0,~] = CalderonCalculusTest(uinc,graduinc,g1p,g1m);
[~,beta1] = CalderonCalculusTest(uinc,graduinc,g2p,g2m);
beta0=-beta0; beta1=-beta1;

% INTEGRAL EQUATIONS AND POTENTIAL REPRESENTATIONS

% System of BIE

lambda1phi2 = [V1(s)       R1*K(s)*R2';...
               R2*J(s)*R1' -W2(s)]\[beta0;beta1];

% Separation of densities by obstacle 

N1 = size(beta0,1); %=length(g1indices)
N2 = size(beta1,1); %=length(g2indices)
lambda1 = lambda1phi2(1:N1);
phi2    = lambda1phi2(N1+1:end);

% Potentials and errors

[SL1,~] = HelmholtzPotentials(g1,xobs);
[~,DL2] = HelmholtzPotentials(g2,xobs);
Uscat = SL1(s)*lambda1+DL2(s)*Q2*phi2;

% ERRORS (accummulate on an existing list if externally run)

ePot = max(abs(Uscat-u(xobs(:,1),xobs(:,2))));

if exist('external')
    errors = [errors; ePot];
else
    errors = [ePot]
end


% PLOT OF THE SOLUTION

if exist('external')
    yesplot=0;
else
    yesplot = input('Do you want a plot? (1 for yes) ');
end

if yesplot
    Box=[-2 5 1.5 -1.5];
    [X,Y,T]=triangulateGeometry(Box,g,0.1,0.1,0.08,0);
    [SL1,~] = HelmholtzPotentials(g1,[X Y]);
    [~,DL2] = HelmholtzPotentials(g2,[X Y]);
    Uscat = SL1(s)*lambda1+DL2(s)*Q2*phi2;
    trisurf(T(:,1:3),X,Y,abs(Uscat))
    view(2), axis equal, axis off, shading interp
    title('Real part of total wave')
end
