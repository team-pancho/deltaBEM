function U=RKCQequation(F,g,k,flag,varargin)

% U=RKCQequation(F,g,k,stages)
% U=RKCQequation(F,g,k,stages,A)
% Input:
%   F(z)    : d x d matrix-valued function of the variable z
%   g       : d*S x (N+1) matrix (S=2 if RadauIIa is used)
%   k       : time step
%  stages   : 1, return all stages, 0 return last stage (step) only
%  varargin : include methods other than 2-stage Radau iia (A) 
% Output:
%   U    : d*S x (N+1) matrix (stages=1)
%   U    : d x (N+1) matrix (stages=0)
%
% Last modified: June 10, 2014

if nargin == 4
    A=[5/12 -1/12; 3/4 1/4];
else
    A=varargin{1};
end

S   = size(A,1);
Am1 = inv(A);
B   = (Am1*ones(size(A,1),1))*[zeros(1,S-1),1]; 

d=size(F(1),1);
N  = size(g,2)-1;
omega= exp(2*pi*1i/(N+1));
R = eps^(0.5/(N+1));

g=bsxfun(@times,g,R.^(0:N));
g=fft(g,[],2);

idx=@(s) (s-1)*d+1:s*d;
    
U=zeros(d*S,N+1);
parfor l=0:floor((N+1)/2) %compute half the hermitian sequence
    [P,Lambda]=eig(Am1-R*omega^(-l)*B);
    Lambda=diag(Lambda)/k;
    gl=kron(inv(P),eye(d))*g(:,l+1);
    ul=zeros(S*d,1);
    for s=1:S
        ul(idx(s))=F(Lambda(s))\gl(idx(s));
    end
    U(:,l+1)=kron(P,eye(d))*ul;
end

%mirror the hermitian sequence
U(:,N+2-(1:floor(N/2)))=conj(U(:,2:floor(N/2)+1));

U=ifft(U,[],2);
U=bsxfun(@times,U,R.^(-(0:N)));
U=real(U);

if flag
    return
end
U=U((S-1)*d+1:end,:);
return
