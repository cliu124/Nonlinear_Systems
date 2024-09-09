clc; clear all; close all;
Nlist2 = []; times = []; rand('state',0);
Nlist = 15;
for N=Nlist
    disp(sprintf('N=%0.0f',N)); rand('state',0);
    metric = @(x_) metric_sphere(x_, 0.2, 0.5, [eps eps 1-eps], 1/N, 1/N/10); 
    %metric = @(x_) metric_circ(x_, 0.1, 0.25, [eps eps], 1/N, 1/N/10); 
    %metric = @(x_) metric_shock(x_, 0.1, 0.1, [0.5 0.5 0.5], 5*1/N, 1/N/10); 
    %metric = @(x_) metric_uniform(x_, [1/N 1/N 1/N]);
    options = gen_options();
    options.qualM = 1;
    %options.qualP = -pi/2*1.2;
    %options.prag_adapt = 1;
    options.debug = 0;
    options.greedyCLR = 2;
    options.smpRFN = 1;
    options.consRFN = 2;
    options.consRM = 0;
    options.innerit = 40;
    options.RMnd3D = 1;
    options.nosparse = false;
    %options.prag_crs = true;
    %options.swap3D = false;
    options.minA = 0;
    rand('state',0);
    %[xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(2*N,options,0); 
    [xy,tri,bndmesh,options,geomfunc] = mesh_sphere(2*N,options,0); 
    tic;
    for i=1:options.outerit
	Nmetric = metric(xy);
	Nmetric = Nmetric*2.^(1./3.);
	[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options); 
    end;
    times(end+1) = toc/size(xy,1)*1e3
    Nlist2(end+1,:) = [size(xy,1) size(tri,1)];
end;
if numel(Nlist) > 1
plot(Nlist/Nlist(end)*size(xy,1)/1e3,times,'.b-'); ylabel('seconds/1k nodes'); xlabel('k nodes')
else
trimesh(bndmesh.fac,xy(:,1),xy(:,2),xy(:,3),'edgecolor','k')
end;

for i=1:9
	options.qualM = i;
	quality = elem_qual(tri,xy,Nmetric,options);
	disp(sprintf('metric %0.0f: average=%0.1f%%, min=%0.1f%%',i,mean(quality)*100,min(quality)*100));
end;
[fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2] = bks_all3D(tri);
options.qualM = 10;
edgQ = elem_qual(edg,xy,Nmetric,options);
facQ = elem_qual(fac,xy,Nmetric,options);
triQ = elem_qual(tri,xy,Nmetric,options);
tmp = quantile(edgQ,[0.25 0.75 0.01 0.99 0 1]); disp(sprintf('edges: 50 %% in [%0.1f %0.1f]%0.1f, 98 %% in [%0.1f %0.1f]%0.1f, all in [%0.1f %0.1f]%0.1f',[tmp(1) tmp(2) tmp(2)/tmp(1) tmp(3) tmp(4) tmp(4)/tmp(3) tmp(5) tmp(6) tmp(6)/tmp(5)]));
tmp = quantile(facQ,[0.25 0.75 0.01 0.99 0 1]); disp(sprintf('areas: 50 %% in [%0.1f %0.1f]%0.1f, 98 %% in [%0.1f %0.1f]%0.1f, all in [%0.1f %0.1f]%0.1f',[tmp(1) tmp(2) tmp(2)/tmp(1) tmp(3) tmp(4) tmp(4)/tmp(3) tmp(5) tmp(6) tmp(6)/tmp(5)]));
tmp = quantile(triQ,[0.25 0.75 0.01 0.99 0 1]); disp(sprintf('vol.s: 50 %% in [%0.1f %0.1f]%0.1f, 98 %% in [%0.1f %0.1f]%0.1f, all in [%0.1f %0.1f]%0.1f',[tmp(1) tmp(2) tmp(2)/tmp(1) tmp(3) tmp(4) tmp(4)/tmp(3) tmp(5) tmp(6) tmp(6)/tmp(5)]));