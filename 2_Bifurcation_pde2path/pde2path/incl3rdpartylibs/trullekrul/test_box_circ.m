clc; clear all; close all;
times = []; rand('state',0);
Nlist = 15; %40:15:100;  %%110:20:190; 
for N=Nlist
    disp(sprintf('N=%0.0f',N));
    metric = @(x_) metric_circ(x_, 0.2, 0.5, [eps eps], 1/N, 1/N/10);
    %metric = @(x_) metric_shock(x_, 0.1, 0.1, [0.5 0.5], 5*1/N, 1/N/10);
    %metric = @(x_) metric_uniform(x_, [1/N 1/N]);
    options = gen_options();
    options.fastRM = 2;
    options.minA = eps;
    options.mntn_bks = 1;
    %options.fastRFN = 0;
    rand('state',0);
    %[xy,tri,bndmesh,options,geomfunc] = mesh_rect_hole(N,0.75,0.1,options);
    %[xy,tri,bndmesh,options,geomfunc] = mesh_circ(2*N,options);
    %[xy,tri,bndmesh,options,geomfunc] = mesh_rect(2*N,options);
    [xy,tri,bndmesh,options,geomfunc] = mesh_rect(2*N,options,1);
    tic;
    for i=1:options.outerit
	Nmetric = metric(xy);
	[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
    end;
    times(end+1) = toc/size(xy,1)*1e3
end;
if numel(Nlist) > 1
     plot(Nlist/Nlist(end)*size(xy,1)/1e3,times,'.b-'); ylabel('seconds/1k nodes'); xlabel('k nodes'); title(sprintf('pragmatic takes %0.2f seconds/1k nodes (for 250k)',24/250)); %1e3/times(end)/(250e3/24) 20% of pragmatic speed
else
	trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','.'); xlim([-0.1 1.1]); ylim([-0.1 1.1]);
end;

for i=1:7
	options.qualM = i;
	quality = elem_qual(tri,xy,Nmetric,options);
	disp(sprintf('metric %0.0f: average=%0.1f%%, min=%0.1f%%',i,mean(quality)*100,min(quality)*100));
end;

[edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg,nd2tri] = bks_all(tri);
options.qualM = 10;
edgQ = elem_qual(edg,xy,Nmetric,options);
triQ_ = elem_qual(tri,xy,Nmetric,options);
tmp = quantile(edgQ,[0.25 0.75 0.01 0.99 0 1]); disp(sprintf('edges: 50 %% in [%0.1f %0.1f]%0.1f, 98 %% in [%0.1f %0.1f]%0.1f, all in [%0.1f %0.1f]%0.0f',[tmp(1) tmp(2) tmp(2)/tmp(1) tmp(3) tmp(4) tmp(4)/tmp(3) tmp(5) tmp(6) tmp(6)/tmp(5)]));
tmp = quantile(triQ_,[0.25 0.75 0.01 0.99 0 1]); disp(sprintf('vol.s: 50 %% in [%0.1f %0.1f]%0.1f, 98 %% in [%0.1f %0.1f]%0.1f, all in [%0.1f %0.1f]%0.1f',[tmp(1) tmp(2) tmp(2)/tmp(1) tmp(3) tmp(4) tmp(4)/tmp(3) tmp(5) tmp(6) tmp(6)/tmp(5)]));