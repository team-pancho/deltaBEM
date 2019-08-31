% Time harmonic elastic scattering by three circle-shaped obstacles.

% Last modified January 12, 2015

%% CHANGE THESE LINES DEPENDING ON WHETHER THE EXAMPLE IS NEW OR NOT

newexample = 0;
nameoffile = '../meshes/ThreeCircles';

%% Simulation Parameters

NP = 100;      % Number of discretization points per obstacle
k  = 3;
s  = -1i*k;    % Wave number
erho = 50;     % Elastic Density 
elambda = 5;   % Lame's first parameter
emu = 3;       % Elastic shear modulus

% Plane pressure elastic wave (incident wave)

ed = [1,0]; ed = ed/(norm(ed));   % Propagation direction
cT = sqrt(emu/erho);              % Shear wave speed
cL = sqrt((elambda+2*emu)/erho);  % Pressure wave speed

u = @(x,y) exp(-(s/cL)*(ed(1)*x+ed(2)*y));
u1 = @(x,y) ed(1)*u(x,y);
u2 = @(x,y) ed(2)*u(x,y);
u1x = @(x,y) -ed(1)^2*s/cL*u(x,y);
u1y = @(x,y) -ed(1)*ed(2)*s/cL*u(x,y);
u2x = @(x,y) -ed(1)*ed(2)*s/cL*u(x,y);
u2y = @(x,y) -ed(2)^2*s/cL*u(x,y);

% Components of the stress tensor
sig11 = @(x,y) 2*emu*u1x(x,y)+elambda*(u1x(x,y)+u2y(x,y));
sig12 = @(x,y) emu*(u1y(x,y)+u2x(x,y));
sig22 = @(x,y) 2*emu*u2y(x,y)+elambda*(u1x(x,y)+u2y(x,y));

%% Geometry

% The obstacles: Three circular obstacles

[g1,g1p,g1m]=sample(@ellipse,NP,[1 1],[1 1]);
[g2,g2p,g2m]=sample(@ellipse,NP,[1 1],[3 3]);
[g3,g3p,g3m]=sample(@ellipse,NP,[1 1],[3.5 .4]);
[g,gp,gm]=merge({g1,g2,g3},{g1p,g2p,g3p},{g1m,g2m,g3m});


%% Meshing

if newexample
    Box=[-1 5.5 5 -1.5];
    h = 0.075;                  % Mesh diameter
    [X,Y,T]=triangulateGeometry(Box,g,0.015,0.015,h,0);
    T=T(:,1:3);
    save(nameoffile,'X','Y','T');
else
    load(nameoffile);
end


%% Calderon Calculus Matrices

% General CC matrices. The flag 1 indicates the use of mixing matrices
[Q,M,~,~] = CalderonCalculusMatrices(g,1);

% Resizing the quadrature matrix (Q) and the mass matrix (M) for elasticity
O = zeros(size(M));
Q = [Q O; O Q];
M = [M O; O M];

% Calderon Calculus Operators. Only J and W are needed
[~,~,J,W] = CalderonCalculusElasticWave(g,gp,gm,emu,elambda,erho);

% Layer potentials
[SL,DL] = ElasticWavePotentials(g,[X,Y],emu,elambda,erho);

%% Boundary Data
% Capturing the measurements of the incident wave on the boundary of the
% scatterers.

[b0x,b1x] = CalderonCalculusTest(...
                          u1,@(x1,x2) [sig11(x1,x2) sig12(x1,x2)],gp,gm,1);
          
[b0y,b1y] = CalderonCalculusTest(...
                          u2,@(x1,x2) [sig12(x1,x2) sig22(x1,x2)],gp,gm,1);
           
beta0 = [b0x;b0y];   % Trace of the displacement
beta1 = [b1x;b1y];   % Normal traction

%% Solving the integral equations for the scattered wave field

% A) Scattered wave. Direct formulation. Dirichlet BC.
        
phiD = -M\beta0;
RHS = -W(s)*phiD;
lambdaD = (.5*M + J(s))\RHS;
UD =  -SL(s)*lambdaD +DL(s)*(Q*phiD);
UD = reshape(UD,.5*size(UD,1),2);

% B) Scattered Wave. Direct formulation. Neumann BC.

lambdaN = -M\beta1;
phiN = -W(s)\((.5*M+J(s))*lambdaN);
UN = -SL(s)*lambdaN+DL(s)*(Q*phiN);
UN = reshape(UN,.5*size(UN,1),2);  

%% Reconstructing the total wave-fields
% The total wave is obtained by superposing the scattered wave to the
% incident wave at the observation points.


Uinc = [u1(X,Y) u2(X,Y)];   % Sampling the incident wave
UDtotal = UD + Uinc;                    % Total Dirichlet wave
UNtotal = UN + Uinc;                    % Total Neumann Wave

%% Plotting
% The norms of the real and imaginary parts of the displacement vector are
% plotted.

realDnorm = sqrt(real(UDtotal(:,1)).^2 + real(UDtotal(:,2)).^2);
imagDnorm = sqrt(imag(UDtotal(:,1)).^2 + imag(UDtotal(:,2)).^2);
realNnorm = sqrt(real(UNtotal(:,1)).^2 + real(UNtotal(:,2)).^2);
imagNnorm = sqrt(imag(UNtotal(:,1)).^2 + imag(UNtotal(:,2)).^2);

subplot(2,2,1)
trisurf(T,X,Y,realDnorm)
title('Dirichlet Problem, Real Part')
shading interp, view(2), colormap copper, axis off

subplot(2,2,2)
trisurf(T,X,Y,imagDnorm)
title('Dirichlet Problem, Imaginary Part')
shading interp, view(2), colormap copper, axis off

subplot(2,2,3)
trisurf(T,X,Y,realNnorm)
title('Neumann Problem, Real Part')
shading interp, view(2), colormap copper, axis off

subplot(2,2,4)
trisurf(T,X,Y,imagNnorm)
title('Neumann Problem, Imaginary Part')
shading interp, view(2), colormap copper, axis off
