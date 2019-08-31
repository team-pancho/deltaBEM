function  g = starshape(N,ep,r,rp)

% g = starshape(N,ep,r,rp)
% Input:
%    N       : discretization parameter
%    ep      : epsilon parameter
%    r       : 2-Pi periodic radius function
%   rp       : first derivative of r
% Output:
%    g       : discrete geometry
%
% Last modified: August 2, 2012

h = 1/N;
t = h*(0:N-1); t = t+ep*h; t=t(:);
tau  = 2*pi*(t-0.5*h);
t    = 2*pi*t;
cost = cos(t);
sint = sin(t);
rt   = r(t);
rpt  = rp(t);

g.midpt = bsxfun(@times,[cost sint], rt);
g.brkpt = bsxfun(@times,[cos(tau) sin(tau)], r(tau));
xp = 2*pi*[rpt.*cost-rt.*sint,...
           rpt.*sint+rt.*cost];
g.normal  = h*[xp(:,2) -xp(:,1)];
g.next    = [2:N 1];
g.comp    = [1];

return