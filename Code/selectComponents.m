function [indices,restr] = selectComponents(g,v)

% [indices,restr] = selectComponents(g,v)
% Input:
%       g   :   joined geometry composed of several geometries
%       v   :   a row vector with numbers corresponding to geometries
% Output:
%       indices  :  a row vector containing the indices corresponding to
%                   the geometries selected
%       restr    : restriction matrix
% Last Modified: August 2 2013

indices = [];
N = length(g.midpt);
g.comp=[g.comp, N+1];
for i = v
    indices = [indices g.comp(i):(g.comp(i+1)-1)];
end
restr=sparse(1:length(indices),indices,1,length(indices),N);
return