% A script to test four formulations for the wave equation
% Last modified: December 4, 2015

% discretization parameters
if ~exist('external')
    N = input('Number of points in space: ');
    M = input('Number of time steps: ');
    T = input('Final time: ');
    Tlag = input('Time lag: ');
    nexp = input('Experiment number: ');
    fork = input('Fork (1 or 0): '); 
end

kappa = T/M;

% BDF2 method transfer function
p =@(z) (1-z)+0.5*(1-z)^2;

% signal
signal  =@(t) (t>=0).*sin(2*t).^5;
signalp =@(t) 10*(t>=0).*(sin(2*t).^4).*cos(2*t);

% direction of wave propogation
dir = [1 1]/sqrt(2);

% Discrete geometry
g  = ellipse(N,0,[1 2],[0 0]);
gp = ellipse(N,1/6,[1 2],[0 0]);
gm = ellipse(N,-1/6,[1 2],[0 0]);

% Generate the signal and boundary data
[Q,MM,Pp,Pm] = CalderonCalculusMatrices(g,fork);
[uincp,dnuincp] = planewave(signal,signalp,dir,Tlag,gp.midpt,gp.normal,T,M,1);
[uincm,dnuincm] = planewave(signal,signalp,dir,Tlag,gm.midpt,gm.normal,T,M,1);
beta0 = (Pp*uincp+Pm*uincm);
beta1 = (Pp*dnuincp+Pm*dnuincm);

% Solve, postprocess, and compute errors
[V,K,J,W] = CalderonCalculusHelmholtz(g,gp,gm,fork);
obs = [0.1 0.2; -0.1 -0.3; 0.5 0];
[S,D] = HelmholtzPotentials(g,obs);
Uex = planewave(signal,signalp,dir,Tlag,obs,0*obs,T,M,0);

switch nexp 
    case  1        
        % Indirect Dirichlet
        eta     = CQequation(V,beta0,kappa,p);
        uIndDir = CQforward(S,eta,kappa,p);
        err     = max(abs(Uex(:,end)-uIndDir(:,end)));
    case  2
        % Indirect Neumann
        psi     = -CQequation(W,beta1,kappa,p);
        uIndNeu = CQforward(D,psi,kappa,p);
        err     = max(abs(Uex(:,end)-uIndNeu(:,end)));
    case  3 
        % Direct Dirichlet
        RHSdir    = 0.5*MM*beta0 + CQforward(K,beta0,kappa,p);
        lambda    = CQequation(V,RHSdir,kappa,p);
        uDirDir   = CQforward(S,lambda,kappa,p) - CQforward(D,Q*beta0,kappa,p);
        errU      = max(abs(Uex(:,end)-uDirDir(:,end)));
        errLambda = N*max(abs(lambda(:,end)-beta1(:,end)));
        err       = [errU errLambda];
    case  4
        % Direct Neumann
        RHSneu  = -0.5*MM'*beta1 + CQforward(J,beta1,kappa,p);
        phi     = -CQequation(W,RHSneu,kappa,p);
        uDirNeu = CQforward(S,beta1,kappa,p)-CQforward(D,Q*phi,kappa,p);
        errU    = max(abs(Uex(:,end)-uDirNeu(:,end)));
        errPhi  = max(abs(phi(:,end)-beta0(:,end)));
        err     = [errU errPhi];
end

if exist('external')
    errors = [errors; err];
else
    errors = err;
end