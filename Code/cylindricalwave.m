function [u,dnu]=cylindricalwave(f,c,x0,x,normal,T,M,option)

% [u,dnu]=cylindricalwave(f,c,[x0 y0],x,normal,T,M,option)
%
% Input:
%      f      : vectorized function of one variable or 
%               cell array of vectorized functions of one variable or
%               Nsrc x M+1 matrix of a sampled signal
%      c      : wave speed
%      [x0 y0]: source point or matrix of source points (corresponding
%               to one source point per handle in f.
%      x      : N x 2 matrix (observation points)
%      normal : N x 2 matrix (vectors for directional observations)
%      T      : final time
%      M      : number of time steps
%      option : 1 (compute everything) 2 (compute only uinc)
% Output:
%      u      : N x (M+1) matrix
%      dnu    : N x (M+1) matrix
%
% Last modified: May 29, 2014

N=size(x,1);
Nsc=size(x0,1);

if isa(f,'function_handle');
    upsample=1;
    k=5;        %upsampling parameter
    t=linspace(0,T,k*M+1);
    F=f(t);
    F=repmat(F,Nsc,1);
    kappaCQ=T/(k*M);
elseif iscell(f)
    upsample=1;
    k=5;        %upsampling parameter
    F=zeros(Nsc,k*M+1);
    t=linspace(0,T,k*M+1);
    for j=1:length(f)
        fj=f{j};
        F(j,:)=fj(t);
    end
    kappaCQ=T/(k*M);
elseif ismatrix(f)
    upsample=0;
    kappaCQ=T/M;
    F=f;
end

u=zeros(N,size(F,2));

for j=1:Nsc
    diffs=bsxfun(@minus,x,x0(j,:));
    dist =sqrt(diffs(:,1).^2+diffs(:,2).^2);
    U=@(s) 1i/4*besselh(0,1,1i*s*dist);
    u=u+CQforward(@(s) U(s/c),F(j,:),kappaCQ);
end

if option==2
    dnu=[];
    if upsample
        u=u(:,1:k:end);
    end
    return
end

dnu=zeros(N,size(F,2));

for j=1:Nsc
    diffs=bsxfun(@minus,x,x0(j,:));
    dist =sqrt(diffs(:,1).^2+diffs(:,2).^2);
    dipoles=sum(diffs.*normal,2)./dist;
    dnU=@(s) s/4*besselh(1,1,1i*s*dist).*dipoles;
    dnu=dnu+CQforward(@(s)dnU(s/c),F(j,:),kappaCQ);
end

if upsample
    u=u(:,1:k:end);
    dnu=dnu(:,1:k:end);
end

return