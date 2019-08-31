function gnew=affine(g,A,b)

% gnew = affine(g,[a11 a12; a21 a22],[b2;b2])
% Input:
%     g    : discrete sampled geometry
%     A    : 2 x 2 matrix
%     b    : 2 x 1 vector
% Output:
%     gnew : discrete sampled geometry
%            corresponding the the affine transformation xnew=A*x+b
% Last modified: August 2, 2013 

if det(A)<=0
    disp('det A <=0 non supported')
    gnew=g;
end

gnew.midpt=bsxfun(@plus,b',g.midpt*A');
gnew.brkpt=bsxfun(@plus,b',g.brkpt*A');
gnew.normal=g.normal*([0 -1; 1 0]*A'*[0 1;-1 0]);
gnew.next=g.next;
gnew.comp=g.comp;
if isfield(g,'parity')
    gnew.parity = g.parity;
end

return

