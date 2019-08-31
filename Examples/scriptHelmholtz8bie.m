% Script to test 8 equations listed in Documentation, including 8 of the
% equations in the paper "A Nystrom flavored Calderon Calculus of order 
% three for two dimensional waves"

% January 12, 2015

if ~exist('external')
    k = input('Input frequency: ');
    N = input('Number of points: ');
    nexp = input('Experiment number (1-8): ');
    fork = input('Use fork? (0 or 1): ');
end

% FREQUENCY

s=-1i*k;

% EXACT SOLUTION

xsc = [0.1 0.2];
r   = @(x,y) sqrt((x-xsc(1)).^2+(y-xsc(2)).^2);
u   = @(x,y) besselh(0,1,1i*s*r(x,y)); % Exterior radiating solution
ux  = @(x,y) -(x-xsc(1)).*1i*s.*besselh(1,1,1i*s*r(x,y))./r(x,y);
uy  = @(x,y) -(y-xsc(2)).*1i*s.*besselh(1,1,1i*s*r(x,y))./r(x,y);
gradu = @(x,y) [ux(x,y) uy(x,y)];

% INCIDENT WAVE

uinc     = @(x,y) -u(x,y);
graduinc = @(x,y) -gradu(x,y);

% DISCRETE GEOMETRY
% tv-shaped domain and ellipse (first semiaxes, then center)

[g1,g1p,g1m]=sample(@tvshape,2*N);
[g2,g2p,g2m]=sample(@ellipse,N,[1 2],[4 5]);
[g,gp,gm]=merge({g1,g2},{g1p,g2p},{g1m,g2m});

% OBSERVATION POINTS
 
xobs = [0 4; 4 0; -4 2; 2 -4];  %(exterior to domains)

% OPERATORS, POTENTIALS & RHSs
    
[Q,M] = CalderonCalculusMatrices(g,fork);
[V,K,J,W] = CalderonCalculusHelmholtz(g,gp,gm,fork);
[beta0,beta1] = CalderonCalculusTest(uinc,graduinc,gp,gm,fork); 
beta0=-beta0; beta1=-beta1;
[SL,DL] = HelmholtzPotentials(g,xobs); % note DL has no Q on right
          
% EXACT SOLUTION BETA_0 AND BETA_1 (only to compute errors)

beta0e = u(g.midpt(:,1),g.midpt(:,2));
beta1e = sum(gradu(g.midpt(:,1),g.midpt(:,2)).*g.normal,2);  
    
% INTEGRAL EQUATIONS AND POTENTIAL REPRESENTATIONS
    
switch nexp
	case 1    % iD01
        eta = V(s)\beta0;
        phi = M\beta0;
        phi = Q*phi;
        lambda = -0.5*eta+M\(J(s)*eta);
        uh = SL(s)*eta;
    case 2    % iD02
        psi = (0.5*M+K(s))\beta0;
        phi = M\beta0;
        phi = Q*phi;
        lambda = -M\(W(s)*psi);
        psi = Q*psi;
        uh = DL(s)*psi;
    case 3    % dD01
        phi = M\beta0;
        lambda = V(s)\(-0.5*M*phi+K(s)*phi);
        phi = Q*phi;
        uh = DL(s)*phi-SL(s)*lambda;
    case 4    % dD02
    	phi = M\beta0;
        lambda = -(0.5*M + J(s))\(W(s)*phi);
        phi = Q*phi;
        uh = DL(s)*phi-SL(s)*lambda;
    case 5    % iN01
        eta  = (-0.5*M + J(s))\beta1;
        phi = M\(V(s)*eta);
        phi = Q*phi;
        lambda = M\beta1;
        uh = SL(s)*eta;
    case 6    % iN02
        psi = - W(s)\beta1;
        phi = 0.5*psi+M\(K(s)*psi);
        phi = Q*phi;
        lambda = M\beta1;
        psi = Q*psi;
        uh  = DL(s)*psi;
    case 7    % dN01
        lambda = M \ beta1;
        phi  = (-0.5*M+ K(s)) \ (V(s)*lambda);
        phi  = Q*phi;
        uh = DL(s)*phi-SL(s)*lambda;   
    case 8    % dN02
        lambda = M\beta1;
        phi = -W(s)\(0.5*M*lambda+J(s)*lambda);
        phi = Q*phi;
        uh  = DL(s)*phi-SL(s)*lambda;
end

% ERRORS (accummulate on an existing list if externally run)

ePot = max(abs(uh-u(xobs(:,1),xobs(:,2))));
eLambda  = max(abs(beta1e-lambda))*N;
ePhi = max(abs((beta0e-phi)));

if exist('external')
    errors = [errors; ePot eLambda ePhi];
else
    errors = [ePot eLambda ePhi];
end

