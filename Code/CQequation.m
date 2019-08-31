function g=CQequation(F,h,k,varargin)

% g=CQequation(F,h,k)
% g=CQequation(F,h,k,p)
%
% Input:
%    F : transfer function            d x d matrix valued
%    h : time values of RHS           d x (N+1) matrix
%    k : time-step
%    p : transfer function of multistep method
% Output:
%    g : F(\partial_k)^{-1} u         d x (N+1) matrix
%
% Last Modified : June 4, 2014

if nargin==3
    p = @(z) 1.5-2*z+0.5*z.^2;     % The default is BDF2 
else
    p=varargin{1};
end

d = size(h,1);             % number of components of problem
N = size(h,2)-1;           % N = number of time-steps
omega = exp(2*pi*1i/(N+1));
R = eps^(0.5/(N+1));

h = bsxfun(@times,h,R.^(0:N));  % scaling
h = fft(h,[],2);                % dft by columns
g = zeros(d,N+1);
parfor l=0:floor((N+1)/2)       %compute half the sequence
    g(:,l+1)=F(p(R.*omega.^(-l))/k)\h(:,l+1);   % \hat w
end

g(:,N+2-(1:floor(N/2)))=conj(g(:,2:floor(N/2)+1));
                                % mirror the hermitian seq
g=real(ifft(g,[],2));                            % w
g=bsxfun(@times,g,R.^(-(0:N)));                  % g

return