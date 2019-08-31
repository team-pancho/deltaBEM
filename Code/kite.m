function  g = kite(N,ep)

% function  g = kite(N,ep)
% Input:
%        N   = number of space intervals
%        ep  = epsilon parameter
% Output:
%         g  : discrete sampled geometry for a kite-shaped domain
% Last modified: August 2, 2013 


h = 1/N;
t = h*(0:N-1); t = t+ep*h;

g.midpt  = [cos(2*pi*t)+cos(4*pi*t);... 
            2*sin(2*pi*t)]';
g.brkpt  = [cos(2*pi*(t-0.5*h))+cos(4*pi*(t-0.5*h));... 
            2*sin(2*pi*(t-0.5*h))]';
xp= [-2*pi*sin(2*pi*t)-4*pi*sin(4*pi*t);...
      4*pi*cos(2*pi*t)]'; 
g.normal  = h*[xp(:,2) -xp(:,1)];
g.next    = [2:N 1];
g.comp    = [1];
return
