function g = joinGeometry(g1,g2)

% function g = joinGeometry(g1,g2)
% Input
%            g1 : 1st geometry
%            g2 : 2nd geometry
% Output
%             g : merged geometry
% Last Modified: August 2, 2013

N=size(g1.midpt,1);
g.midpt = [g1.midpt; g2.midpt];
g.brkpt = [g1.brkpt; g2.brkpt];
g.normal = [g1.normal; g2.normal];
g.next = [g1.next,N+g2.next];
g.comp = [g1.comp,N+g2.comp];

if (~isfield(g1,'parity'))&&(~isfield(g2,'parity'))
    return
elseif (isfield(g1,'parity'))&&(isfield(g2,'parity'))
    g.parity = blkdiag(g1.parity,g2.parity);
elseif ~isfield(g1,'parity')
    g1.parity = sparse(size(g1.midpt,1),size(g1.midpt,1));
else
    g2.parity = sparse(size(g2.midpt,1),size(g2.midpt,1));
end
g.parity = blkdiag(g1.parity,g2.parity);

return