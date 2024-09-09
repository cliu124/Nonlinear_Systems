clc; clear all; close all;
N = 1000;
xy = rand(N,2);
tri = delaunay(xy(:,1),xy(:,2));
%trimesh(tri,xy(:,1),xy(:,2),'k','linewidth',2); box off; axis off;

flname = 'mesh.pvd';
sclrvr = 1+sin(xy(:,1)*10).*cos(xy(:,2)*10);
ii = 1;
gen_vtk(flname,tri,xy,sclrvr,ii);
sclrvr = 1+sin(xy(:,1)*10).*cos(xy(:,2)*10);
ii = 2;
sclrvr = 1+sin(0.2+xy(:,1)*10).*cos(0.2+xy(:,2)*10);
gen_vtk(flname,tri,xy,sclrvr,ii);