function [xy,tri,bndmesh,options,geomfunc] = mesh_rect_hole(N,dy,r,options)
x = [0 0 r 0 0 1 1 r r]'; y = [0 dy-r dy dy+r 1 1 0 dy+r dy-r]'; 
tri = delaunay(x,y);
tri = elem_fix(tri,[x y]);
xy = [x y];
triy = mean(y(tri),2);
tri = tri(abs(triy-0.75)>1e-6,:);

%options = gen_options();
bndmesh = bndmesh_polygon(tri,xy,[],options);
%generate geomfunc
distC    = @(x_) geom_circ(x_,r,[0 dy]); 
distR    = @(x_) geom_rect(x_); 
geomfunc = @(x_) geom_diff(distR(x_),distC(x_));
%uniform metric
Nmetric = [ones(size(xy,1),1) zeros(size(xy,1),1) ones(size(xy,1),1)]*N;
%adapt
[tri,xy,Nmetric,bndmesh] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
%shuffle nodes
[C,I] = sort(rand(size(xy,1),1));
xy = xy(I,:);
[C,I] = sort(I);
tri = I(tri);
bndmesh.edg = sort(I(bndmesh.edg),2);
options.area = 0;