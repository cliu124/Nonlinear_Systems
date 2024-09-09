clc; clear all; close all;
rand('state',0);
do3D = 1;
N=25;
    disp(sprintf('N=%0.0f',N));
    times = []; %1000.8 1513.6 1022.59    861.37
    metric = @(x_) metric_iso(x_, 1/N); 
    options = gen_options();
    options.prag_adapt = 2;
    options.debug = 0;
    options.qualM = 1;
    options.consRM = 0;
    options.fastRM = 1;
    options.fstRFN = 0;
    options.consRFN = 0;
    options.greedyCLR = 2;
    options.verbose = 1;
    options.mntn_bks = 1;
    options.RMnd3D = 1;
    options.minsliver = false;

    %options.RMedg = 1; options.nosparse = false;
    rand('state',0);
    if do3D
    	options.innerit=10;
    	%options.forpres = '3a';
        %options.qualM = 8;
        %options.minqual = 1e-13;
    	[xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(N,options,2);
    else
    options.forpres = 'a';
    [xy,tri,bndmesh,options,geomfunc] = mesh_rect(N,options,1);
    end;
	Nmetric = metric(xy);
	[tri,xy,Nmetric,bndmesh] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);

    
    if do3D
    	%options.forpres = '3b';
        options.innerit = 30;
    	metric = @(x_) metric_shock(x_, 0.1, -0.1, [0.5 0.5 0.5], 5*1/N, 1/N/10); 
    else
    	options.forpres = 'b';
    	options.innerit=20;
	metric = @(x_) metric_shock(x_, 0.1, -0.1, [0.5 0.5], 5*1/N, 1/N/10); 
    end;
	Nmetric = metric(xy); 
	[tri,xy,Nmetric,bndmesh] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
   
    
% clc; clear all;
% for i=[(0:40) 121]
%     load(sprintf('./3a/%0.0f.mat',i));
%     [fac,fac2tri,tri2fac] = bks_all3D(tri);
%     trimesh(fac(fac2tri(:,1)~=0,:),xy(:,1),xy(:,2),xy(:,3),'edgecolor','k'); axis off; box off;
%     title(intit,'fontsize',32); print('-dpng',sprintf('./3a/%0.0f.png',i));
% end;
% for i=0:121
%     load(sprintf('./3b/%0.0f.mat',i));
%      [fac,fac2tri,tri2fac] = bks_all3D(tri);
%     trimesh(fac(fac2tri(:,1)~=0,:),xy(:,1),xy(:,2),xy(:,3),'edgecolor','k'); axis off; box off;
%     title(intit,'fontsize',32); print('-dpng',sprintf('./3b/%0.0f.png',i));
% end;
    
