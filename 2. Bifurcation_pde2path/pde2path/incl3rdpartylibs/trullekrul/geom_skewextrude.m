function dist = geom_skewextrude(x_,func2D,zrange,skewvec) %we do not take edges at ends into account, so edge function is needed.
x_(:,1:2) = x_(:,1:2) - repmat(skewvec,size(x_,1),1).* repmat((x_(:,3)-zrange(1))/(zrange(2)-zrange(1)),1,2);
dist2D = func2D(x_(:,1:2));
dist = [dist2D zeros(size(x_,1),1)];
ends = or(abs(x_(:,3)-skewvec(1)) < 1e3*eps,abs(x_(:,3)-skewvec(2)) < 1e3*eps);
dist(ends,:) = 0.;


%clc; clear all;
%N = 10; options = gen_options();
%options.consRM = 0;
%options.smpRFN = 1;
%options.consRFN = 2;
%%options.minA = 0;
%options.innerit = 20;
%options.prag_crs = true;
%[xy,tri,bndmesh,options,geomfunc] = mesh_cyl(N,options,0);
%xy(:,1) = xy(:,1) + xy(:,3);
%geomcirc = @(x_) geom_circ(x_,0.5,[0.5 0.5]);
%geomfunc = @(x_) geom_skewextrude(x_,geomcirc,[0. 1.],[1. 0.]);
%geomcirc = @(x_) geom_circ3D(x_,0.5,[0.5 0.5 0.]);
%geomcrc2 = @(x_) geom_move(x_,geomcirc,[1. 0. 1.]);
%geomcrcs = @(x_) geom_union(geomcirc(x_),geomcrc2(x_));
%geomfunc = {geomfunc, geomcrcs};
%Nmetric = metric_uniform(xy, [1/N 1/N 1/N]);
%[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
%
%close all; trimesh(bndmesh.fac,xy(:,1),xy(:,2),xy(:,3));
%tmp = geomfunc{1}(xy(bndmesh.fac(:),:)); max(abs(tmp(:,1)))
%[crnds,edg1,edg2] = geom_crnds(bndmesh);
%tmp = geomfunc{2}(xy(edg1(:),:)); max(abs(tmp(:,1)))