function [XX,YY,TT,LBN] = separateTriangulation(X,Y,T,sub,gBn)

% [XX,YY,TT,LBN] = separateTriangulation(X,Y,T,sub,gBn)
%
% Input:
%
%   X,Y   : nVert x 1 vectors with coordinates of vertices
%   T     : nElt x 4 with triangulation (4th column is subdomain).
%   sub   : Vector containing the ordering of the subdomains.
%   gBn   : Vector containing the -global- list of boundary nodes.(optional)
%  
% Output:
%
%     XX, YY : cell array with X,Y coordinates separated by subdomain.
%     TT     : cell array with triangulation separated by subdomain.
%     LBN    : cell array containing the -local- list of boundary nodes
%              separated by subdomain. (optional only available if gBn is
%              provided as input).
%
% Last modified: August 29, 2014.

nVert = size(X,1);
nSubd = max(T(:,4));

for j = 1:nSubd
    
    i = find(T(:,4)==sub(j));
    TT{j} = T(i,1:3); 
    vert = TT{j}(:);
    vert = unique(vert);
    XX{j} = X(vert);
    YY{j} = Y(vert);
    transp = zeros(nVert,1); 
    transp(vert) = 1:length(vert);
    TT{j} = transp(TT{j});
    
    % Optional list of boundary nodes
    if nargin == 5 
        LBN{j} = vert(ismember(vert,gBn));
        LBN{j} = transp(LBN{j});
    end  
    
end

end

    
    
    