function U=RKCQforward(F,g,k,flag,varargin)

% U=RKCQforward(F,g,k,flag)
% U=RKCQforward(F,g,k,flag,A)
% Input:
%   F(z)    : d2 x d1 matrix-valued function of the variable z
%   g       : S*d1 x (N+1) matrix (S=2 if RadauIIa is used)
%   k       : time step
%  flag   : 1, return all stages, 0 return last stage (step) only
%  varargin : include methods other than 2-stage Radau iia (A)    
% Output:
%   U    : S*d2 x (N+1) matrix (stages=1)
%   U    : d2 x (N+1) matrix (stages=0)
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

[d1,d2] = size(F(1));
N  = size(g,2)-1;
omega= exp(2*pi*1i/(N+1));
R = eps^(0.5/(N+1));

g=bsxfun(@times,g,R.^(0:N));
g=fft(g,[],2);

idx=@(s) (s-1)*d1+1:(s*d1);
idy=@(s) (s-1)*d2+1:(s*d2); 

U=zeros(S*d1,N+1);

parfor l=0:floor((N+1)/2) %compute half the hermitian sequence
    [P,Lambda]=eig(Am1-R*omega^(-l)*B);
    Lambda=diag(Lambda)/k;
    gl=kron(inv(P),eye(d2))*g(:,l+1);
    ul=zeros(S*d1,1);
    for s=1:S
        ul(idx(s))=F(Lambda(s))*gl(idy(s));
    end
    U(:,l+1)=kron(P,eye(d1))*ul;
end

%mirror the hermitian sequence
U(:,N+2-(1:floor(N/2)))=conj(U(:,2:floor(N/2)+1));

U=ifft(U,[],2);
U=bsxfun(@times,U,R.^(-(0:N)));
U=real(U);

if flag
    return
end
U=U((S-1)*d1+1:end,:);
return


