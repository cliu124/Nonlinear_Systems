function faceplot(p,u) 
% faceplot: plot u on faces of 3D domain
%
% faceplot(p,u)
%
try; p.pdeo.grid.plotFaces(u,'EdgeColor',p.plot.EdgeColor); catch; p.pdeo.grid.plotFaces(u); end 
set(gca,'FontSize',p.plot.fs); colormap(p.plot.cm);


