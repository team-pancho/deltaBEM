function [X,Y,T,sub,gBn] = triangulateGeometry(r, g,d_abs,d_rel,hgrid,varargin)

% [X,Y,T,sub,gBn]=triangulateGeometry([xmin xmax ymax ymin],g,dabs,drel,hgrid,fill)
%
% Input:
%
% [xmin xmax ymax ymin] : framing rectangle
% g                     : geometry of the curve (discrete geometry)
% dabs                  : absolute distance of the grid to the curve
% drel                  : relative distance of the grid to the curve
% hgrid                 : meshgrid size
% fill                  : 1 generate grid for the interior
%                         0 or missing -> interior not meshed
%
% Output:
%
%   X,Y   : Npt x 1 vectors with X and Y coordinates of points
%   T     : Nelements x 4 matrix with elements (nodes and subdomains)
%   sub   : Subdomain ordering 
%   gBn   : Column vector containig the -global- list of boundary nodes.
%           
%
% Note that you need the PDE toolbox in your Matlab distribution
%
% Last modified: August 29, 2014.

if ~isempty(varargin)
    intQ=varargin{1};
    intQ=double(intQ);
    intQ=intQ(1);
else
    intQ=0;
end

[pde_fig,ax]=pdeinit;

n_curves=length(g.comp);
n=length(g.midpt);
index=[g.comp n+1];
d_abs=abs(d_abs); d_rel=abs(d_rel);

lengths=sqrt(g.normal(:,1).^2+g.normal(:,2).^2);

counter=1;
pderect(r,'R1');
Op='R1'; Op2=[];

for j=1:n_curves
    counter=counter+1;
    label2=['R' num2str(counter)];
    indAux=[index(j):index(j+1)-1];
    curve=g.midpt(indAux,:)+...
        bsxfun(@times,g.normal(indAux,:),d_abs./lengths(indAux)+d_rel);
    pdepoly(curve(:,1).', curve(:,2).', label2);
    Op=[Op '-' label2];
    if intQ==1
        curve=g.midpt(indAux,:)-...
            bsxfun(@times,g.normal(indAux,:),d_abs./lengths(indAux)+d_rel);
        counter=counter+1;
        label2=['R' num2str(counter)];
        pdepoly(curve(:,1).', curve(:,2).',label2);
        Op2=[Op2 '+' label2]; 
        object(j,:)=curve(1,1:2);
    end
end


% Generating the mesh
Op=['(' Op ')' Op2];

pdetool('appl_cb',1);
set(ax,'XLim',[r(1) r(2)]);
set(ax,'YLim',[r(4) r(3)]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String',Op)

% Mesh generation:
setappdata(pde_fig,'trisize',hgrid);
setappdata(pde_fig,'Hgrad',1.25+hgrid);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
pdetool('initmesh')

% Extracting the information about the triangulation
h = findobj(get(pde_fig,'Children'),'flat','Tag','PDEMeshMenu');

hp = findobj(get(h,'Children'),'flat','Tag','PDEInitMesh');
p = get(hp,'UserData');  % Node coordinates

ht = findobj(get(h,'Children'),'flat','Tag','PDEMeshParam');
t = get(ht,'UserData');  % Delauney triangulation

he=findobj(get(h,'Children'),'flat','Tag','PDERefine');
e = get(he,'UserData');   % Nodes that define the (boundary) edges
gBn = unique(e(1:2,:));   % Boundary Nodes


X=p(1,:)';
Y=p(2,:)';
T=t';


if intQ==1
    for j=1:n_curves
        vertex=intersect(find(object(j,1)==X),find(object(j,2)==Y));
        [elt,v]=find(T(:,1:3)==vertex);
        sub(j)=T(elt(1),4);
    end
    ext = setdiff(1:n_curves+1,sub);
    sub=[ext sub];
else
    sub=1;
end

return
