function frequencyDomainPlot(X,Y,Tri,U,Nsnap,name,g,yes)

% frequencyDomainPlot(X,Y,Tri,U,Nsnap,name,g,yes)
% Input:
%     [X,Y,Tri] : triangulation in the plane
%      U        : complex values of a function on (X,Y)
%      Nsnap    : number of desired snapshots
%      name     : name of file for plot
%      g        : scatterer
%      yes      : 1 (include the scatterers), 0 (do not)
%
% Last modified: August 6, 2013

% Names of files

indices={'01','02','03','04','05','06','07','08','09',...
         '10','11','12','13','14','15','16','17','18','19',...
         '20','21','22','23','24','25','26','27','28','29',...
         '30','31','32','33','34','35','36','37','38','39',...
         '40','41','42','43','44','45','46','47','48','49',...
         '50','51','52','53','54','55','56','57','58','59',...
         '60','61','62','63','64','65','66','67','68','69',...
         '70'};
for i=1:length(indices)
    longname{i}=[name,indices{i},'.png'];
end
Nsnap=min(Nsnap,length(indices));

% Unpacking the geometries and fixing heights

G=unpackGeometry(g);
nComp=length(G);

Umax=max(real(U)); Umin=min(real(U));
U(1:2)=[];
T=linspace(0,2*pi,Nsnap+1); T(end)=[];

% Running times

i=0;
for t=T
    trisurf(Tri,X,Y,...
            [Umax;Umin;real(U)*cos(t)+imag(U)*sin(t)]);
    view(2), shading interp, axis equal, axis off
    hold on
    if yes
        for k=1:nComp
            plot3(G{k}.midpt(:,1),G{k}.midpt(:,2),0*G{k}.midpt(:,1),'-');
        end
    end
    i=i+1;
    saveas(gcf,longname{i})
    hold off
    pause(0.02)
end
    

