% Script to test 2 boundary integral equations iD01 and iN02 for the 
% exterior Laplace equation. 

% January 12, 2015

if ~exist('external')
    N = input('Number of points: ');
    nexp = input('Experiment number 1 or 2: ');
    fork = input('Use fork? (0 or 1): ');
end

% EXACT SOLUTION

u  = @(x,y) x./(x.^2+y.^2);
ux = @(x,y) 1./(x.^2+y.^2)-2*x.^2./(x.^2+y.^2).^2;
uy = @(x,y) -2*x.*y./(x.^2+y.^2).^2;
gradu = @(x,y) [ux(x,y) uy(x,y)];

% DISCRETE GEOMETRY
% tv-shaped domain and ellipse (first semiaxes, then center)

[g1,g1p,g1m]=sample(@tvshape,2*N);
[g2,g2p,g2m]=sample(@ellipse,N,[1 2],[4 5]);
[g,gp,gm]=merge({g1,g2},{g1p,g2p},{g1m,g2m});

% OBSERVATION POINTS
 
xobs = [0 4; 4 0; -4 2; 2 -4];  %(exterior to domains)

% OPERATORS, POTENTIALS, RHSs

[Q,M] = CalderonCalculusMatrices(g,fork);
[V,K,J,W,C] = CalderonCalculusLaplace(g,gp,gm,fork);
[beta0,beta1] = CalderonCalculusTest(u,gradu,gp,gm,fork); 
[SL,DL] = LaplacePotentials(g,xobs);

% EXACT BETA_0 AND BETA_1

beta0e = u(g.midpt(:,1),g.midpt(:,2));
beta1e = sum(gradu(g.midpt(:,1),g.midpt(:,2)).*g.normal,2);  

% INTEGRAL REPRESENTATION AND POTENTIAL REPRESENTATIONS

one = ones(size(g.midpt(:,1)));
zero = zeros(size(g.midpt(:,1)));

switch nexp
	case 1    % iD01
        eta = [ V one ; one' 0 ]\[ beta0; 0 ];
        eta(end) = [];
        phi = M\beta0;
        phi = Q*phi;
        lambda = -0.5*eta+M\(J*eta);
        uh = SL*eta;
    case 2    % iN02
        psi = - (W + C)\beta1;
        phi = 0.5*psi+M\(K*psi);
        phi = Q*phi;
        lambda = M\beta1;
        psi = Q*psi;
        uh  = DL*psi;
end

% ERRORS (accumulate on an existing list if externally run)

ePot = max(abs(uh-u(xobs(:,1),xobs(:,2))));
eLambda  = max(abs(beta1e-lambda))*N;
ePhi = max(abs((beta0e-phi)));

if exist('external')
    errors = [errors; ePot eLambda ePhi];
else
    errors = [ePot eLambda ePhi];
end
