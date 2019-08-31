function  g = ellipse(N,ep,R,c)

% g = ellipse(N,eps,[a,b],[cx,cy])
% Input:
%    N       : discrete interval number
%    ep      : epsilon parameter
%    [a b]   : semiaxes
%    [cx,cy] : center
% Output:
%    g       : sampled discrete geometry
%
% Last modified: August 2, 2013

h = 1/N;
t = h*(0:N-1); t = t+ep*h;
cost = cos(2*pi*t);
sint = sin(2*pi*t);

g.midpt = [c(1)+R(1)*cost;...
           c(2)+R(2)*sint]';
g.brkpt = [c(1)+R(1)*cos(2*pi*(t-0.5*h));...
           c(2)+R(2)*sin(2*pi*(t-0.5*h))]';         
xp = [-R(1)*2*pi*sint;...
       R(2)*2*pi*cost]';
g.normal  = h*[xp(:,2) -xp(:,1)];
g.next    = [2:N 1];
g.comp    = [1];

return