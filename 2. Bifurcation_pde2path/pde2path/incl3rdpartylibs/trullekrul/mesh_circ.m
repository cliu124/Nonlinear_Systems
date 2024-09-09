function [xy,tri,bndmesh,options,geomfunc] = mesh_circ(N,options,ustruct)
if nargin ~= 2
[xy,tri,bndmesh,options,geomfunc] = mesh_rect(N,options,ustruct);
else
[xy,tri,bndmesh,options,geomfunc] = mesh_rect(N,options);
end; x = xy(:,1); y = xy(:,2);
%transform cube to sphere
bndmesh = bndmesh_polygon(tri,[x y],[],options);
x = 2*(x - 0.5); y = 2*(y - 0.5);
xy = [x.*sqrt(1.- 0.5*y.^2) ...
      y.*sqrt(1.- 0.5*x.^2)];
xy = xy/2 + 0.5;
bndmesh.crnds = []; bndmesh.IDs(:) = 1;
options.area = 0;
geomfunc = @(x_) geom_circ(x_,0.5,[0.5 0.5]);