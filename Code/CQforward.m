function u = CQforward(F,g,k,varargin)

% u = CQforward(F,g,k)
% u = CQforward(F,g,k,p)
%
% Input:
%    F : transfer function (d1 x d2 matrix-valued)
%    g : time values of input (d2 x (N+1) matrix)
%    k : time-step
%    p : transfer function of multistep method
% Output:
%    u : F(\partial_k) g              d1 x (N+1) matrix
%
% Last Modified : June 4, 2014

if nargin==3
    p = @(z) 1.5-2*z+0.5*z.^2;     % The default is BDF2 
else
    p=varargin{1};
end

d1 = size(F(1),1);             % number of components of output vector
N  = size(g,2)-1;              % N = number of time-steps

omega = exp(2*pi*1i/(N+1));
R = eps^(0.5/(N+1));

h = bsxfun(@times,g,R.^(0:N));    
h = fft(h,[],2);               % DFT by columns (\hat H)
u = zeros(d1,N+1);
parfor l=0:floor((N+1)/2)
    u(:,l+1)=F(p(R*omega^(-l))/k)*h(:,l+1);     % \hat v   
end
u(:,N+2-(1:floor(N/2)))=conj(u(:,2:floor(N/2)+1));

u=real(ifft(u,[],2));                           % v
u=bsxfun(@times,u,R.^(-(0:N)));

return