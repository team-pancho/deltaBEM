function  g = tvshape(N,ep)

% function  g = tvshape(N,ep)
% Input:
%        N   : number of space intervals
%        ep  : epsilon parameter
% Output:
%        g   : discrete sampling of smoothened square
%
% Last modified: August 2, 2013 

h = 1/N;
t = h*(0:N-1); t = t+ep*h;

g.midpt  = [(1+cos(2*pi*t).^2).*cos(2*pi*t);...
            (1+sin(2*pi*t).^2).*sin(2*pi*t)]';
g.brkpt  = [(1+cos(2*pi*(t-0.5*h)).^2).*cos(2*pi*(t-0.5*h));...
            (1+sin(2*pi*(t-0.5*h)).^2).*sin(2*pi*(t-0.5*h))]';
R=[cos(pi/4) sin(pi/4); -sin(pi/4) cos(pi/4)]';
g.midpt = g.midpt*R;
g.brkpt = g.brkpt*R;

xp = [6.*pi.*sin(2.*pi.*t).^3 - 8.*pi.*sin(2.*pi.*t);...
      2.*pi.*(4.*cos(2.*pi.*t) - 3.*cos(2.*pi.*t).^3)]';
xp = xp*R;
g.normal  = h*[xp(:,2) -xp(:,1)];
g.next    = [2:N 1];
g.comp    = [1];
return
