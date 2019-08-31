%a script to thoroughly test RKCQ methods with various interior problems

%last modified March 29, 2016

if ~exist('external')
    N = input('Number of points: ');
    M = input('Number of time steps: ');
    T = input('Final Time: ');
    nexp = input('Experiment number (1-2): ');
    tlag = input('Time lag: ');
end

% Set up problem parameters
kappa=T/M;

%Radau iia method
A=[5/12 -1/12; 3/4 1/4];
b=A(end,:);
St=size(A,1);
c=A*ones(St,1);

%geometry
g=kite(N,0);
gp=kite(N,1/6);
gm=kite(N,-1/6);

%signal data
signal=@(t) sin(2*t).^9.*(t>=0);
signalp=@(t) 18*sin(2*t).^8.*cos(2*t).*(t>=0);
dir=[-1 1]/sqrt(2);


%observation points (interior)
obs=[0 -1.25; 1.5 0];

%incident wave observations
[Q,Mass,~,~] = CalderonCalculusMatrices(g);
QKron=kron(eye(St),Q);

beta0=zeros(N,St,M);
beta1=zeros(N,St,M);

for st=1:St
    [up,dnup]=planewave(@(t) signal(t+kappa*c(st)),@(t) signalp(t+kappa*c(st)),...
        dir,tlag,gp.midpt,gp.normal,T-kappa,M-1,1);
    [um,dnum]=planewave(@(t) signal(t+kappa*c(st)),@(t) signalp(t+kappa*c(st)),...
        dir,tlag,gm.midpt,gm.normal,T-kappa,M-1,1);
    beta0(:,st,:) = 0.5*up+0.5*um;
    beta1(:,st,:) = Q*(0.5*dnup+0.5*dnum);
end

beta0=reshape(beta0,[St*N,M]);
beta1=reshape(beta1,[St*N,M]);

% Calderon calculus
[V,K,J,W] = CalderonCalculusHelmholtz(g,gp,gm);
MassKron=kron(eye(St),Mass);
[S,D]=HelmholtzPotentials(g,obs);

switch nexp
        
    case 1

        % indirect Neumann Problem with RKCQequation
        tic;
        psi=RKCQequation(W,-beta1,kappa,1,A);
        Uh=RKCQforward(D,QKron*psi,kappa,0,A);
        Uh=[zeros(size(obs,1),1) Uh];
        Uex = planewave(signal,signalp,dir,tlag,obs,0*obs,T,M,0);
        ePot=max(abs(Uh(:)-Uex(:)));
        toc
        
    case 2
        
        % indirect Dirichlet Problem with RKCQequation
        tic;
        eta=RKCQequation(V,beta0,kappa,1,A);
        Uh=RKCQforward(S,eta,kappa,0,A);
        Uh=[zeros(size(obs,1),1) Uh];
        Uex = planewave(signal,signalp,dir,tlag,obs,0*obs,T,M,0);
        ePot=max(abs(Uh(:)-Uex(:)));
        toc
        
end

% ERRORS (accummulate on an existing list if externally run)


if exist('external')
    error = [error; ePot];
else
    error = ePot
end