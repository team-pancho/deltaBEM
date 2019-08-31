function [u,v,sxx,sxy,syy]=sourceSolutionLame(mu,lambda,x0,d)

% [u,v,sxx,sxy,syy]=sourceSolutionLame(mu,lambda,[x0,y0],[d1 d2])
% Input:
%     mu,lambda : Lame parameters
%     [x0,y0]   : source point
%     [d1,d2]   : direction of displacement vector
% Output:
%    Five vectorized functions of (x,y), corresponding to a source solution
%    of the Lame equation and its corresponding stress tensor
% Last modified: October 31, 2013

xx = x0(1);
yy = x0(2);
d1 = d(1);
d2 = d(2);

rr = @(x,y) (x-xx).^2+(y-yy).^2;
rd = @(x,y) (x-xx)*d1+(y-yy)*d2;

A  = mu/(lambda+2*mu);
C1 = -(1+A)/(4*pi*mu);
C2 = (1-A)/(4*mu*pi);

u = @(x,y) C1*0.5*log(rr(x,y))*d1...
           +C2*rd(x,y)./rr(x,y).*(x-xx);
v = @(x,y) C1*0.5*log(rr(x,y))*d2...
           +C2*rd(x,y)./rr(x,y).*(y-yy);
ux = @(x,y) (C1+C2)*d1*(x-xx)./rr(x,y)...
            -2*C2*rd(x,y)./(rr(x,y).^2).*(x-xx).^2 ...
            +C2*rd(x,y)./rr(x,y);
uy = @(x,y) C1*d1*(y-yy)./rr(x,y)...
            +C2*d2*(x-xx)./rr(x,y)...
            -2*C2*rd(x,y)./(rr(x,y).^2).*(x-xx).*(y-yy);
vx = @(x,y) C1*d2*(x-xx)./rr(x,y)...
            +C2*d1*(y-yy)./rr(x,y)...
            -2*C2*rd(x,y)./(rr(x,y).^2).*(y-yy).*(x-xx);
vy = @(x,y) (C1+C2)*d2*(y-yy)./rr(x,y)...
            -2*C2*rd(x,y)./(rr(x,y).^2).*(y-yy).^2 ...
            +C2*rd(x,y)./rr(x,y);

sxx = @(x,y) 2*mu*ux(x,y)+lambda*(ux(x,y)+vy(x,y));
sxy = @(x,y) mu*(uy(x,y)+vx(x,y));
syy = @(x,y) 2*mu*vy(x,y)+lambda*(ux(x,y)+vy(x,y));
return
       