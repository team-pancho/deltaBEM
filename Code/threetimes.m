function [gnew,gpnew,gmnew]=threetimes(utility,g,gp,gm,varargin)

% [gnew,gpnew,gmnew]=threetimes(@utility,g,gp,gm,varargin)
% Input:
%     @utility : handle to one of the geometric utilities
%     varargin : other parameters needed by utility
% Output:
%     gnew     : utility(g,varargin)
%     gpnew    : utility(gp,varargin)
%     gmnew    : utility(gm,varargin)
% Last modified: January 12, 2015

switch nargin
    case 4
        gnew =utility(g);
        gpnew=utility(gp);
        gmnew=utility(gm);
    case 5
        gnew =utility(g,varargin{1});
        gpnew=utility(gp,varargin{1});
        gmnew=utility(gm,varargin{1});
    case 6
        gnew =utility(g,varargin{1},varargin{2});
        gpnew=utility(gp,varargin{1},varargin{2});
        gmnew=utility(gm,varargin{1},varargin{2});
end
end