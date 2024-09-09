clc; clear all; close all;
Nlist2 = []; times = []; rand('state',0);
Nlist = 5;
for N=Nlist
    disp(sprintf('N=%0.0f',N)); rand('state',0);
    metric = @(x_) metric_sphere(x_, 0.2, 0.5, [eps eps 1-eps], 1/N, 1/N/10); 
    %metric = @(x_) metric_circ(x_, 0.1, 0.25, [eps eps], 1/N, 1/N/10); 
    %metric = @(x_) metric_shock(x_, 0.1, 0.1, [0.5 0.5 0.5], 5*1/N, 1/N/10); 
    metric = @(x_) metric_uniform(x_, [1/N 1/N 1/N]);
    options = gen_options();
    options.qualM = 1;
    %options.qualP = -pi/2*1.2;
    %options.prag_adapt = 1;
    options.debug = 0;
    options.greedyCLR = 2;
    options.smpRFN = 1;
    options.consRFN = 2;
    options.consRM = 0;
    options.verbose = 0;
    options.innerit = 20;
    options.RMnd3D = 1;
    options.nosparse = false;
    options.prag_crs = true;
    options.swap3D = false;
    options.verbose = 1;
    options.mntn_bks = 1;
    options.minA = 0;
%    options.fastRFN = false;
%    options.nosparse = false;
%    options.minA = 1e-15;
    rand('state',0);
    %[xy,tri,bndmesh,options,geomfunc] = mesh_rect_hole3D(N,[12 4 4],1,options);
    [xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(2*N,options,0); 
    %[xy,tri,bndmesh,options,geomfunc] = mesh_sphere(2*N,options,0); 
    %bndmesh = bndmesh_polyhedron(tri,xy,bndmesh,options); bndmesh.IDs(and(and(xy(bndmesh.fac(:,1),1) == 0,xy(bndmesh.fac(:,2),1) == 0),xy(bndmesh.fac(:,3),1) == 0)) = -1;
    tic;
    for i=1:options.outerit
	Nmetric = metric(xy); %try
	Nmetric = Nmetric*2.^(1./3.);
	[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options); 
	%end;
    end;
    times(end+1) = toc/size(xy,1)*1e3
    Nlist2(end+1,:) = [size(xy,1) size(tri,1)];
end;
if numel(Nlist) > 1
%trixy = (xy(tri(:,1),:)+xy(tri(:,2),:)+xy(tri(:,3),:)+xy(tri(:,4),:))/4; I = sqrt(sum((trixy-repmat([0.5 0.5 0.5],size(tri,1),1)).^2,2))<0.25; min(triQ(I))
%save('uniform.mat','Nlist2','times');
%save('sphere.mat','Nlist2','times');
plot(Nlist/Nlist(end)*size(xy,1)/1e3,times,'.b-'); ylabel('seconds/1k nodes'); xlabel('k nodes') %size(xy,1)/times(end)/(250e3/24) 10% of pragmatic speed
else	
[fac,fac2tri,tri2fac] = bks_all3D(tri); figure(); facxy = reshape(xy(fac(:),:),size(fac,1),9); facxy = [mean(facxy(:,1:3),2) mean(facxy(:,4:6),2) mean(facxy(:,7:9),2)]; minmax = [min(xy); max(xy)];
I = abs(facxy(:,1)-minmax(1,1)) < 1e-15; subplot(2,3,1); trimesh(fac(I,:),xy(:,2),xy(:,3)); title('xmin'); xlim(minmax(:,2)'+[-0.1 0.1]); ylim(minmax(:,3)'+[-0.1 0.1]); 
I = abs(facxy(:,1)-minmax(2,1)) < 1e-15; subplot(2,3,4); trimesh(fac(I,:),xy(:,2),xy(:,3)); title('xmax'); xlim(minmax(:,2)'+[-0.1 0.1]); ylim(minmax(:,3)'+[-0.1 0.1]); 
I = abs(facxy(:,2)-minmax(1,2)) < 1e-15; subplot(2,3,2); trimesh(fac(I,:),xy(:,1),xy(:,3)); title('ymin'); xlim(minmax(:,1)'+[-0.1 0.1]); ylim(minmax(:,3)'+[-0.1 0.1]); 
I = abs(facxy(:,2)-minmax(2,2)) < 1e-15; subplot(2,3,5); trimesh(fac(I,:),xy(:,1),xy(:,3)); title('ymax'); xlim(minmax(:,1)'+[-0.1 0.1]); ylim(minmax(:,3)'+[-0.1 0.1]); 
I = abs(facxy(:,3)-minmax(1,3)) < 1e-15; subplot(2,3,3); trimesh(fac(I,:),xy(:,1),xy(:,2)); title('zmin'); xlim(minmax(:,1)'+[-0.1 0.1]); ylim(minmax(:,2)'+[-0.1 0.1]); 
I = abs(facxy(:,3)-minmax(2,3)) < 1e-15; subplot(2,3,6); trimesh(fac(I,:),xy(:,1),xy(:,2)); title('zmax'); xlim(minmax(:,1)'+[-0.1 0.1]); ylim(minmax(:,2)'+[-0.1 0.1]); 
%	trimesh(tri,xy(:,1),xy(:,2),xy(:,3));

end;

for i=[(1:7) 8 9]
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


%xy = [0 0 sqrt(2/3); ... 
      %-sqrt(1/3) 0 0; ...
      %cos(pi/3)*sqrt(1/3) sin(pi/3)* sqrt(1/3) 0; ...
      %cos(pi/3)*sqrt(1/3) sin(-pi/3)*sqrt(1/3) 0];
%Nmetric = repmat([1 0 1 0 0 1],4,1);
%xy = [-sqrt(1/3) 0; ...
      %cos(pi/3)*sqrt(1/3) sin(pi/3)* sqrt(1/3); ...
      %cos(pi/3)*sqrt(1/3) sin(-pi/3)*sqrt(1/3)];
%Nmetric = repmat([1 0 1],3,1);
%
%
%tri = repmat(1:size(xy,2)+1,2,1); 
%L1 = sqrt(sum((xy(tri(:,1),:)-xy(tri(:,2),:)).^2,2))
%L2 = sqrt(sum((xy(tri(:,3),:)-xy(tri(:,2),:)).^2,2))