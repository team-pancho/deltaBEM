function [g,gp,gm]=merge(G,Gp,Gm)

% [g,gp,gm]=merge({g1,g2,...,gM},{gp1,gp2,...,gpM},{gm1,...,gmM});
% Input:
%     g1,g2,...   : geometry data structures (collected in cell array)
%     gp1,gp2,... : (same)
%     gm1,gm2,... : (same)
% Output:
%     g  : geometry data structure 
%          join(join(....(join(g1,g2),g3),...),gM) [joinGeometry]
% Last modified: January 12, 2015

M=length(G);
g =G{1};
gp=Gp{1};
gm=Gm{1};
for c=2:M
    g =joinGeometry(g,G{c});
    gp=joinGeometry(gp,Gp{c});
    gm=joinGeometry(gm,Gm{c});
end
end

    