function [V,K,J,W] = CalderonCalculusElasticWave(g,gp,gm,mu,lambda,rho)

% [V,K,J,W] = CalderonCalculusElasticWave(g,gp,gm,mu,lambda,rho)
% Input:
%     g : principal geometry
%    gp : companion geometry with epsilon = 1/6
%    gm : companion geometry with epsilon = -1/6
%    mu,lambda : Lame parameters
%    rho : density
% Output:
%     V : single layer operator (function of s) 
%     K : double layer operator K 
%     J : transpose double layer operator K 
%     W : hypersingular operator W
%
% Last Modified: March 29, 2016

[Vp,Kp,Jp,WRp,WSp] = CalderonCalculusEWHalf(g,gp,mu,lambda,rho);
[Vm,Km,Jm,WRm,WSm] = CalderonCalculusEWHalf(g,gm,mu,lambda,rho);

[Q,~,Pp,Pm] = CalderonCalculusMatrices(g,1);  % Fork = 1
O = zeros(size(Q));
Q  = [Q O; O Q]; 
Pp = [Pp O;O Pp];
Pm = [Pm O;O Pm];

V = @(s) Pp*Vp(s)+Pm*Vm(s);
K = @(s) (Pp*Kp(s)+Pm*Km(s))*Q;
J = @(s) Q*(Pp*Jp(s)+Pm*Jm(s));
W = @(s) Q*(Pp*WRp(s)+Pm*WRm(s))*Q + Pp*WSp(s)+Pm*WSm(s);

end

% Subfunction computing the two halves of the operators

function [V,K,J,WR,WS] = CalderonCalculusEWHalf(g,gp,mu,lambda,rho)

% [V,K,J,WR,WS] = CalderonCalculusEWHalf(g,gp,mu,lambda, rho)
% Input:
%     g : principal geometry
%    gp : companion geometry
%    mu,lambda : Lame parameters
%   rho : density
% Output:
%     V : single layer operator (function of s)
%     K : double layer operator
%     J : transpose double layer operator
%     WR: regular part of the hypersingular operator 
%     WS: singular part of the hypersingular operator
%
% Last Modified: July 15, 2014.

cT = sqrt(mu/rho);
cL = sqrt((lambda+2*mu)/rho);
nu = .5*lambda/(lambda+mu); 
CC = 2*(1-nu);  % CC=(lambda+2*mu)/(lambda+mu)
xi = sqrt(mu/(lambda+2*mu));

% Some basic blocks

R1  = bsxfun(@minus,gp.midpt(:,1),g.midpt(:,1)');
R2  = bsxfun(@minus,gp.midpt(:,2),g.midpt(:,2)');
R = sqrt(R1.^2+R2.^2); 
R11 = R1.*R1;
R12 = R1.*R2;
R21 = R12;
R22 = R2.*R2;
RaRb = [R11 R12; R21 R22];
O = zeros(size(R));

% Utilities

Id = @(A) [A O; O A];
Sc = @(A) [A A; A A];

% Four functions

psi  = @(r) besselk(0,r/cT)...
          +(cT./r).*(besselk(1,r/cT)-xi*besselk(1,r/cL));
psi_r = @(r) -(1/cT)*besselk(1,r/cT)...
             -(2*cT./r.^2).*(besselk(1,r/cT)-xi*besselk(1,r/cL))... 
             -(1./r).*(besselk(0,r/cT)-xi^2*besselk(0,r/cL));   
         
chi   = @(r) besselk(2,r/cT)-xi^2*besselk(2,r/cL);
chi_r = @(r) -0.5/cT*(besselk(1,r/cT)+besselk(3,r/cT)...
                        -xi^3*(besselk(1,r/cL)+besselk(3,r/cL)));
                    
% Single layer operator

V = @(s) 1/(2*pi*mu)*(Id(psi(s*R))-Sc(chi(s*R)./R.^2).*RaRb);

% Double layer operator

R1N1 = bsxfun(@times,R1,g.normal(:,1)');
R1N2 = bsxfun(@times,R1,g.normal(:,2)');
R2N1 = bsxfun(@times,R2,g.normal(:,1)');
R2N2 = bsxfun(@times,R2,g.normal(:,2)');
RN   = R1N1+R2N2;
RaNb = [R1N1 R1N2; R2N1 R2N2];
NaRb = [R1N1 R2N1; R1N2 R2N2];

      
K = @(s) -s/(2*pi)*Sc(psi_r(s*R)./R).*(Id(RN)+NaRb+(lambda/mu)*RaNb)...
          +1/(2*pi)*Sc(chi(s*R)).*...
                (-4*Sc(RN./R.^4).*RaRb + Sc(1./R.^2).*NaRb...
                 +Id(RN./R.^2) + (2+lambda/mu)*Sc(1./R.^2).*RaNb)...
          +s/(2*pi)*Sc(chi_r(s*R)).*...
                (2*Sc(RN./R.^3).*RaRb+lambda/mu*Sc(1./R).*RaNb);

% Transposed double layer operator

R1 = -R1; R2 = -R2;  % For J,  r = y-x 
RbRa= RaRb;

N1R1 = bsxfun(@times,gp.normal(:,1),R1);
N2R1 = bsxfun(@times,gp.normal(:,2),R1);
N1R2 = bsxfun(@times,gp.normal(:,1),R2);
N2R2 = bsxfun(@times,gp.normal(:,2),R2);
NR   = N1R1+N2R2;
NbRa = [N1R1 N2R1; N1R2 N2R2]; 
RbNa = [N1R1 N1R2; N2R1 N2R2];

J = @(s) -s/(2*pi)*Sc(psi_r(s*R)./R).*(Id(NR)+NbRa+(lambda/mu)*RbNa)...
          +1/(2*pi)*Sc(chi(s*R)).*...
                (-4*Sc(NR./R.^4).*RbRa+Sc(1./R.^2).*NbRa...
                 +Id(NR./R.^2)+(2+lambda/mu)*Sc(1./R.^2).*RbNa)...
          +s/(2*pi)*Sc(chi_r(s*R)).*...
                (2*Sc(NR./R.^3).*RbRa+lambda/mu*Sc(1./R).*RbNa);
             
% Hypersingular operator: auxiliary computations

N11 = gp.normal(:,1)*g.normal(:,1)';
N12 = gp.normal(:,1)*g.normal(:,2)';
N21 = gp.normal(:,2)*g.normal(:,1)';
N22 = gp.normal(:,2)*g.normal(:,2)';
NaNb = [N11 N12; N21 N22];
NbNa = [N11 N21; N12 N22];
NdotN = N11+N22;

M1 = [N11.*R11+N12.*R21+R11.*N11+R12.*N21,...
            N11.*R12+N12.*R22+R11.*N12+R12.*N22;...
        N21.*R11+N22.*R21+R21.*N11+R22.*N21,...
            N21.*R12+N22.*R22+R21.*N12+R22.*N22];
M2 = [N11.*R11+N12.*R21+R11.*N11+R12.*N21,...
            N21.*R11+N22.*R21+R21.*N11+R22.*N21;...
        N11.*R12+N12.*R22+R11.*N12+R12.*N22,...
            N21.*R12+N22.*R22+R21.*N12+R22.*N22];
M3   = N11.*R11+N12.*R12+N21.*R21+N22.*R22;    

Gp = @(r) -1/(2*pi*rho*cT)*(besselk(1,r/cT) - xi*besselk(1,r/cL));
       
Gpp = @(r) 1/(4*pi*rho*cT^2)*...
             ( besselk(0,r/cT)+besselk(2,r/cT)...
                -xi^2*(besselk(0,r/cL)+besselk(2,r/cL)) );
          
Gppp = @(r) -1/(8*pi*rho*cT^3)*...
            ( 3*besselk(1,r/cT)+besselk(3,r/cT)...
               -xi^3*(3*besselk(1,r/cL)+besselk(3,r/cL)) );
           
Gpppp = @(r) 1/(2*pi*rho*cT^4)*...
               ( (3*cT^2./r.^2+1).*besselk(2,r/cT)...
                  -xi^4*(3*cL^2./r.^2+1).*besselk(2,r/cL));
LapG   = @(r) Gp(r)./r+Gpp(r);
BiLapG = @(r) Gpppp(r)+2*Gppp(r)./r-Gpp(r)./r.^2+Gp(r)./r.^3;
a = @(r) Gpp(r)-Gp(r)./r;
b = @(r) Gp(r)./r;
              
% Regular part of the hypersingular operator

HessF = @(s) Sc(a(s*R)./R.^2).*RaRb +Id(b(s*R));
WR = @(s) CC*s^2*...
           (mu*Sc(BiLapG(s*R)).*(lambda*NaNb + mu*(NbNa + Id(NdotN)))...                  
            -(1/cL^2)*(...
                lambda^2*Sc(LapG(s*R)).*NaNb...
                +2*lambda*mu*(Sc(a(s*R)./R.^2).*M1+2*Sc(b(s*R)).*NaNb)...
                +mu^2*(Id(a(s*R)./R.^2.*M3+b(s*R).*NdotN)...
                       +Sc(NdotN).*HessF(s)...
                       +Sc(a(s*R)./R.^2).*M2+2*Sc(b(s*R)).*NbNa)));
 
% Principal (singular) part of the hypersingular operator
 
R1  = bsxfun(@minus,gp.brkpt(:,1),g.brkpt(:,1)');
R2  = bsxfun(@minus,gp.brkpt(:,2),g.brkpt(:,2)');
R   = sqrt(R1.^2+R2.^2); 

HessFbrk = @(s) Sc(a(s*R)./R.^2).*[R1.*R1 R1.*R2; R2.*R1 R2.*R2]...
                 +Id(b(s*R));
              
Nelt = size(g.brkpt,1);
Dt = -speye(Nelt)+sparse(1:Nelt,g.next,1);
Dt = [Dt O; O Dt];
WS = @(s) Dt*( 4*mu^2*(Id(LapG(s*R)) - HessFbrk(s)) )*Dt';  
      
end        



