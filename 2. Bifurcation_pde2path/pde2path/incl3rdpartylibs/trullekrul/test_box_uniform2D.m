%clc; clear all; close all; N = 15
times = []; rand('state',0);
metric = @(x_) metric_uniform(x_, [1/N 1/N]);
options = gen_options();
options.fastRM = 2;
options.verbose = 0;
rand('state',0);
[xy,tri,bndmesh,options,geomfunc] = mesh_rect(2*N,options,1);
Nmetric = metric(xy);
[tri,xy,Nmetric,~,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,[],options);