function [g,gp,gm]=sample(curve,N,varargin)

% [g,gp,gm]=sample(@curve,N,varargin)
% Input:
%     @curve   : handle to one of the geometry functions
%     N        : number of points
%     varargin : other parameters needed by curve
% Output:
%     g,gp,gm  : sampled geometries with eps=0,1/2,-1/6
% Last modified: January 12, 2015

switch nargin
    case 2
        g =curve(N,0);
        gp=curve(N,1/6);
        gm=curve(N,-1/6);
    case 4
        g =curve(N,0,varargin{1},varargin{2});
        gp=curve(N,1/6,varargin{1},varargin{2});
        gm=curve(N,-1/6,varargin{1},varargin{2});
end
end

