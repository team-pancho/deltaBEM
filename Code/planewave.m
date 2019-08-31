function [u,dnu]=planewave(f,fp,d,tlag,x,normal,T,M,option)

% [u,dnu]=planewave(f,fp,d,tlag,x,normal,T,M,all)
%
% Input:
%      f      : vectorized function of one variable
%      fp     : its derivative
%      d      : 1 x 2 vector (direction of the wave)
%      tlag   : scalar 
%      x      : N x 2 matrix (observation points)
%      normal : N x 2 matrix (vectors for directional observations)
%      T      : final time
%      M      : number of time steps
%      all    : 1 (compute everything) 0 (compute only uinc)
% Output:
%      u      : Nobs x (M+1) matrix
%      dnu    : Nobs x (M+1) matrix
% Last modified: April 23, 2015

t=linspace(0,T,M+1);
dir=sum(bsxfun(@times,x,d),2);        % x . d
u=f(bsxfun(@minus,t-tlag,dir));

if option==0
    dnu=[];
    return
end

nor=sum(bsxfun(@times,normal,d),2);   % normal . d
dnu=fp(bsxfun(@minus,t-tlag,dir));
dnu=-bsxfun(@times,dnu,nor);

return


