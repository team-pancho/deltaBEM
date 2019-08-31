%
% SCRIPT FOR SCATTERING BY (TWO) PENETRABLE OBSTACLES
%
% Last Modified: January 12, 2015

%% CHANGE THESE LINES DEPENDING ON WHETHER THE EXAMPLE IS NEW OR NOT

newexample = 0;
nameoffile = '../meshes/TPtwoScat';

%% Simulation

% Physical parameters and incident wave

c=[1/2 2];
alpha=[1 1];
k=15; s=-1i*k;
dir=[1/sqrt(2) -1/sqrt(2)];   % direction of incident wave

uinc  = @(x,y) exp(s*(dir(1)*x+dir(2)*y));
uincx = @(x,y) s*dir(1)*uinc(x,y); 
uincy = @(x,y) s*dir(2)*uinc(x,y);
graduinc = @(x,y) [uincx(x,y) uincy(x,y)];

% Geometry

N = 600;
N1=N; N2=N; NN=N1+N2;

g1  = tvshape(N1,0);
g1p = tvshape(N1,1/6);
g1m = tvshape(N1,-1/6);

g2  = ellipse(N2,0,[1 1],[3 0]);
g2p = ellipse(N2,1/6,[1 1],[3 0]);
g2m = ellipse(N2,-1/6,[1 1],[3 0]);

g  = joinGeometry(g1,g2);
gp = joinGeometry(g1p,g2p);
gm = joinGeometry(g1m,g2m);

[g1,g1p,g1m]=sample(@tvshape,N1);
[g2,g2p,g2m]=sample(@ellipse,N2,[1 1],[3 0]);
[g,gp,gm]=merge({g1,g2},{g1p,g2p},{g1m,g2m});

if newexample
    gg1  = tvshape(200,0);
    gg2  = ellipse(200,0,[1 1],[3 0]);
    gg  = joinGeometry(gg1,gg2);  % less refined geometry to create triang
end

% From here on, the code ignores how many obstacles there are. All
% information will be used automatically on the geometry, and interior
% domains will be decoupled and separated as needed

% Calderon Calculus

[V,K,J,W]=CalderonCalculusHelmholtz(g,gp,gm);
[Q,M]  =CalderonCalculusMatrices(g);
[Vint,Kint,Jint,Wint]=...
        CalderonCalculusHelmholtzDecoupled(g,gp,gm,c,alpha);
[testu,testdnu] = CalderonCalculusTest(uinc,graduinc,gp,gm);
beta0=-testu;
beta1=-testdnu;

% Costabel-Stephan formulation

tau0=-M\beta0;
tau1=-M\beta1;
Matrix=[W(s)+Wint(s) J(s)+Jint(s);...
       -K(s)-Kint(s) V(s)+Vint(s)];
RHS=[W(s)  0.5*M+J(s); 0.5*M-K(s)  V(s)]*[tau0; tau1];
Sol=Matrix\RHS;
phim   =Sol(1:NN);
lambdam=Sol(NN+1:end);
phip   =-tau0+phim;
lambdap=-tau1+lambdam;

% Postprocessing and separation of components of the solution

phim = Q*phim;
phip = Q*phip;

G = unpackGeometry(g);
for comp=1:length(G)
    N(comp)=size(G{comp}.midpt,1);
end
N=[0 cumsum(N)];
for comp=1:length(G);
    phi{comp}=phim(N(comp)+1:N(comp+1));
    lambda{comp}=(1/alpha(comp))*lambdam(N(comp)+1:N(comp+1));
end


%% Plots

% Triangulations

if newexample
    [X,Y,T,sub]=triangulateGeometry([-2.5 5 2.5 -2.5],gg,0.05,0.05,0.04,1);
    save(nameoffile,'X','Y','T','sub');
else
    load(nameoffile);
end
[XX,YY,TT]=separateTriangulation(X,Y,T,sub);  % Separation of interior/exterior

% Plot of the real part of the solution

Uinc  = uinc(XX{1},YY{1});
[Sext,Dext]=HelmholtzPotentials(g,[XX{1} YY{1}]);
Uscat = Dext(s)*phip-Sext(s)*lambdap;

for sb=1:length(XX)-1
    [S,D]=HelmholtzPotentials(G{sb},[XX{sb+1} YY{sb+1}]);
    Vh{sb} = S(s/c(sb))*lambda{sb}-D(s/c(sb))*phi{sb};
end

trisurf(TT{1},XX{1},YY{1},real(Uinc+Uscat))
view(2), axis equal, axis off, shading interp
hold on
for sb=2:length(XX)
    trisurf(TT{sb},XX{sb},YY{sb},real(Vh{sb-1})), shading interp;
end
saveas(gcf,'figures/HTPrealpartk15.png')

%% Time-harmonic movie (needs to run the Plot part in advance)

Times=linspace(0,2*pi,21); Times(end)=[];  % only one period

i=0;
for t=Times
    trisurf(TT{1},XX{1},YY{1},...
            cos(t)*real(Uinc+Uscat)+sin(t)*imag(Uinc+Uscat))
    view(2), axis equal, axis off, shading interp
    hold on
    for sb=2:length(XX)
        trisurf(TT{sb},XX{sb},YY{sb},...
                cos(t)*real(Vh{sb-1})+sin(t)*imag(Vh{sb-1}))
        shading interp
    end
    hold off
    pause(0.05)
    i=i+1;
    longname=['movies/transmission',num2str(i,'%02d'),'.png'];
    saveas(gcf,longname)
end








    
    
