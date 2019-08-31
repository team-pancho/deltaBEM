function g = openarc(N,ep,x,xp)

% function g = openarc(N,ep,x,x')
% Input:
%       N      :      number of segments
%       ep     :      epsilon parameter
%       x      :      parameterization vector of the curve
%       xp     :      first derivative of x
% Output:
%       g      :      discrete geometry
% Last Modified: August 2, 2013

h = 1/(2*N);
t = h*(0:2*N-1);
t = t+(ep+0.5)*h;

T = @(tau) 0.5*cos(pi*(2*tau-1))+0.5;
TP = @(tau) -pi*sin(pi*(2*tau-1));
xT = x(T(t'));
xpT = xp(T(t'));

g.midpt = xT;
g.brkpt = x(T(t'-h/2));
g.normal = h*[xpT(:,2).*TP(t'), -xpT(:,1).*TP(t')];
g.next = [2:2*N 1];
g.comp = [1];

g.parity = speye(2*N);
g.parity = g.parity-g.parity(end:-1:1,:);
    
return