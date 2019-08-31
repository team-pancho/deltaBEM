function gnew = mirror(g)

% gnew = mirror(g)
% 
% Input : 
%
%    g  : deltaBEM geometry
%
% Output: 
%
%    gnew : geometry mirrored across the x-axis
%
% Last Modified: May 13, 2014

g.brkpt=g.brkpt(g.next,:);

N=size(g.midpt,1);

gnew.midpt=[g.midpt(:,1),-g.midpt(:,2)];
gnew.brkpt=[g.brkpt(:,1),-g.brkpt(:,2)];
gnew.comp=g.comp;
gnew.normal(:,1)=g.normal(:,1);
gnew.normal(:,2)=-g.normal(:,2);

prev=zeros(1,N);

prev(g.next)=1:N;
gnew.next=prev;

return