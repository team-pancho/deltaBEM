function [SL,DL]=ElasticWavePotentials(g,z,mu,lambda,rho)

% [SL,DL]=ElasticWavePotentials(g,z,mu,lambda,rho)
% 
% Input:
%     g : geometry
%     z : K x 2 matrix with points where potentials are evaluated
%     mu, lambda : Lame parameters
%     rho : density
% Output:
%   SL  : matrix valued function of the variable s
%   DL  : matrix valued function of the variable s
%
% Last modified: July 15, 2014.

% Parameters

cT = sqrt(mu/rho);
cL = sqrt((lambda+2*mu)/rho);
xi = sqrt(mu/(lambda+2*mu));    % (cT/cL)


% Basic Matrix blocks

R1  = bsxfun(@minus,z(:,1),g.midpt(:,1)');
R2  = bsxfun(@minus,z(:,2),g.midpt(:,2)');
RaRb = [R1.*R1 R1.*R2; R2.*R1 R2.*R2];
R1N1 = bsxfun(@times,R1,g.normal(:,1)');
R1N2 = bsxfun(@times,R1,g.normal(:,2)');
R2N1 = bsxfun(@times,R2,g.normal(:,1)');
R2N2 = bsxfun(@times,R2,g.normal(:,2)');
RbNa = [R1N1 R2N1; R1N2 R2N2];
RaNb = [R1N1 R1N2; R2N1 R2N2];
RN   = R1N1+R2N2;
R = sqrt(R1.^2+R2.^2);  
O = zeros(size(R));

% Utilities

Id = @(A) [A O; O A]; % Block Identity
Sc = @(A) [A A; A A]; % "Scalar" Matrix

% Auxiliary functions

psi  = @(r) besselk(0,r/cT)...
          +(cT./r).*(besselk(1,r/cT)-xi*besselk(1,r/cL));
      
psi_r = @(r) -(1/cT)*besselk(1,r/cT)...
             -(2*cT./r.^2).*(besselk(1,r/cT)-xi*besselk(1,r/cL))...
             -(1./r).*(besselk(0,r/cT)-xi^2*besselk(0,r/cL));               

chi   = @(r) besselk(2,r/cT)-xi^2*besselk(2,r/cL);

chi_r = @(r) -0.5/cT*...
            (besselk(1,r/cT)+besselk(3,r/cT)...
             -xi^3*(besselk(1,r/cL)+besselk(3,r/cL)));

% Elastic wave Potentials

SL = @(s) 1/(2*pi*mu)*(Id(psi(s*R))-Sc(chi(s*R)./R.^2).*RaRb);

DL = @(s) -s/(2*pi)*Sc(psi_r(s*R)./R).*( Id(RN)+RbNa+(lambda/mu)*RaNb )...
          +1/(2*pi)*Sc(chi(s*R)).*( -4*Sc(RN./R.^4).*RaRb...
                                    +Sc(1./R.^2).*RbNa + Id(RN./R.^2)...
                                    +(2+lambda/mu)*Sc(1./R.^2).*RaNb )...
          +s/(2*pi)*Sc(chi_r(s*R)).*( 2*Sc(RN./R.^3).*RaRb ...
                                       +lambda/mu*Sc(1./R).*RaNb );
             

return