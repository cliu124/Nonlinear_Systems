function [xy,tri,bndmesh,options,geomfunc] = mesh_cyl(N,options,ustruct)
[xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(2*N,options,ustruct);
 x = xy(:,1); y = xy(:,2); z = xy(:,3);
%transform cube to sphere
bndmesh = bndmesh_polyhedron(tri,xy,[],options);
x = 2*(x - 0.5); y = 2*(y - 0.5);
xy = [x.*sqrt(1.- 0.5*y.^2) ...
      y.*sqrt(1.- 0.5*x.^2)];
xy = xy/2 + 0.5;
xy = [xy z];
bndmesh.crnds = []; bndmesh.IDs(:) = 1;
options.area = 0;
cyl      = @(x_) geom_cyl(x_,0.5,[0.5 0.5 0.]);
rect     = @(x_) geom_rect(x_,[1. 1. 1.]);
geomfunc = @(x_) geom_intersect(rect(x_),cyl(x_));
geomcirc = @(x_) geom_circ3D(x_,0.5,[0.5 0.5 0.]);
geomcrc2 = @(x_) geom_move(x_,geomcirc,[0. 0. 1.]);
geomcrcs = @(x_) geom_union(geomcirc(x_),geomcrc2(x_));
geomfunc = {geomfunc, geomcrcs};
bndmesh.IDs(0.5-abs(mean(reshape(xy(bndmesh.fac,3),size(bndmesh.fac)),2)-0.5) < 1e3*eps) = 7;
Nmetric = metric_uniform(xy, [1/N 1/N 1/N]);
[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);

%clc; clear all;
%N = 10; options = gen_options(); ustruct = 0;
%options.consRM = 0;
%options.smpRFN = 1;
%options.consRFN = 2;
%%options.minA = 0;
%options.innerit = 40;
%%options.prag_crs = true;
%close all; trimesh(bndmesh.fac,xy(:,1),xy(:,2),xy(:,3));
%tmp = geomfunc{1}(xy(bndmesh.fac(:),:)); max(abs(tmp(:,1)))
%[crnds,edg1,edg2] = geom_crnds(bndmesh);
%tmp = geomfunc{2}(xy(edg1(:),:)); max(abs(tmp(:,1)))