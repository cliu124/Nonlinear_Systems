%clc; clear all; close all; N = 30; Lx = 1.5; r = 0.01; options = gen_options(); % options.smpRFN = 1; options.consRFN = 2; options.innerit = 30;
function [xy,tri,bndmesh,geomfunc] = mesh_Lfillet(N,Lx,r,options)
x = [0 0 1 1 1+r 1+Lx 1+Lx]'; y = [0 1+Lx 1+Lx 1+r 1 1 0]'; 
tri = [1 2 3; 1 3 4; 1 4 5; 1 5 6; 1 6 7];
tri = elem_fix(tri,[x y]);
xy = [x y];

%options = gen_options();
bndmesh = bndmesh_polygon(tri,xy,[],options);
Nmetric = [ones(size(xy,1),1) zeros(size(xy,1),1) ones(size(xy,1),1)]*1e2;
metric = @(x_) [ones(size(xy,1),1) zeros(size(xy,1),1) ones(size(xy,1),1)]*1e2;
%generate geomfunc
thigh    = @(x_) geom_rect(x_,[1,1+Lx]);
calf     = @(x_) geom_rect(x_,[Lx,1],[1 0]);
L        = @(x_) geom_union(thigh(x_),calf(x_));
cr1      = @(x_) geom_rect(x_,[r,r],[1 1]);
Lr       = @(x_) geom_union(L(x_),cr1(x_));
cr2      = @(x_) geom_circ(x_,r,[1 1]+r); 
geomfunc = @(x_) geom_diff(Lr(x_),cr2(x_));
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
%trimesh(tri,xy(:,1),xy(:,2),'-k');
