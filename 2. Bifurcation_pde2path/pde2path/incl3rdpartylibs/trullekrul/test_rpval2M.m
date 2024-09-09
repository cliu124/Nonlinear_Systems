clc; clear all; N = 30;
%rand('state',0);
[xy,tri,geomfunc] = gen_rect(N); bndmesh = gen_polygon_bndmesh(tri,xy);
nds = tri(:);
tris = [1:size(tri,1) 1:size(tri,1) 1:size(tri,1)]'; 
tic; nd2tri = rpval2M(nds,tris); 
%nd2tri
toc
max(tri(:))
tri(nd2tri(end,1),:)