function L = latticeGeometry(g,centers)

% function L = latticeGeometry(g,centers)
% Input:
%       g       :     the geometry from which the lattice will be made
%       centers :     N by 2 matrix of desired centers for the geometries
% Output:
%       L       :     a single structure of the combined geometries
% Last Modified: July 11, 2013

N = length(centers(:,1));
Larray = [];
for i = 1:N
    Larray = [Larray affine(g,eye(2),centers(i,:)')];
end
l = length(Larray);
L = Larray(1);
for i = 2:l
    L = joinGeometry(L,Larray(i));
end
return
