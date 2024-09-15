clc; %clear all; close all;
times = []; rand('state',0);
Nlist = 2; do3D = true; % Nlist=10; do3D=0; 
for N=Nlist
    disp(sprintf('N=%0.0f',N));
    options = gen_options();
    options.fastRM = 2;
    options.minA = eps;
    options.mntn_bks = 1;
    options.outerit = 5;
    if do3D
        options.smpRFN = 1;
        options.consRFN = 2;
        options.consRM = 0;
        options.outerit = 1;
        [xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(2*N,options,0); 
    else
        %[xy,tri,bndmesh,options,geomfunc] = mesh_rect_hole(N,0.75,0.1,options);
        %[xy,tri,bndmesh,options,geomfunc] = mesh_circ(2*N,options);
        %[xy,tri,bndmesh,options,geomfunc] = mesh_rect(2*N,options);
        [xy,tri,bndmesh,options,geomfunc] = mesh_rect(2*N,options,1);
    end;
    rand('state',0);

    tic;
    for i=1:options.outerit
    r = (sqrt(sum(xy.^2,2))-0.5)/0.02; z = 1./(1+exp(-r)); z=z+10; 
    Nmetric = metric_pnorm(tri,xy,z,0.002); %*0.02^(1-do3D/3)
    %Nmetric, pause
    gen_vtk('pnorm.pvd',tri,xy,z,i);
	[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,[Nmetric z],bndmesh,geomfunc,options);
          %z = Nmetric(:,7); Nmetric = Nmetric(:,1:6); gen_vtk('pnorm.pvd',tri,xy,z,i);
    end;
    r = (sqrt(sum(xy.^2,2))-0.5)/0.02; z = 1./(1+exp(-r));
    z = Nmetric(:,4); Nmetric = Nmetric(:,1:3);
    gen_vtk('pnorm.pvd',tri,xy,z,i+1);
    times(end+1) = toc/size(xy,1)*1e3
end;
if numel(Nlist) > 1
     plot(Nlist/Nlist(end)*size(xy,1)/1e3,times,'.b-'); ylabel('seconds/1k nodes'); xlabel('k nodes'); title(sprintf('pragmatic takes %0.2f seconds/1k nodes (for 250k)',24/250)); %1e3/times(end)/(250e3/24) 20% of pragmatic speed
else
    if do3D
        [fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2] = bks_all3D(tri);
        figure(1); trimesh(bndmesh.fac,xy(:,1),xy(:,2),xy(:,3),'edgecolor','k');
        figure(2); trimesh(bndmesh.fac,xy(:,1),xy(:,2),xy(:,3),z,'facecolor','none'); 
        colormap(winter); axis off; box off; colorbar(); set(gca,'fontsize',24); print('-dpng','pnorm.png'); print('-depsc2','pnorm.eps');
    else
        [edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg,nd2tri] = bks_all(tri);
        figure(1); trisurf(tri,xy(:,1),xy(:,2),z); xlim([-0.1 1.1]); ylim([-0.1 1.1]);
        figure(2); trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','.'); xlim([-0.1 1.1]); ylim([-0.1 1.1]);
    end;
end;

for i=1:7+2*do3D
	options.qualM = i;
	quality = elem_qual(tri,xy,Nmetric,options);
	disp(sprintf('metric %0.0f: average=%0.1f%%, min=%0.1f%%',i,mean(quality)*100,min(quality)*100));
end;

options.qualM = 10;
edgQ = elem_qual(edg,xy,Nmetric,options);
triQ_ = elem_qual(tri,xy,Nmetric,options);
tmp = quantile(edgQ,[0.25 0.75 0.01 0.99 0 1]); disp(sprintf('edges: 50 %% in [%0.1f %0.1f]%0.1f, 98 %% in [%0.1f %0.1f]%0.1f, all in [%0.1f %0.1f]%0.0f',[tmp(1) tmp(2) tmp(2)/tmp(1) tmp(3) tmp(4) tmp(4)/tmp(3) tmp(5) tmp(6) tmp(6)/tmp(5)]));
if do3D
    facQ = elem_qual(fac,xy,Nmetric,options);
    tmp = quantile(facQ,[0.25 0.75 0.01 0.99 0 1]); disp(sprintf('areas: 50 %% in [%0.1f %0.1f]%0.1f, 98 %% in [%0.1f %0.1f]%0.1f, all in [%0.1f %0.1f]%0.1f',[tmp(1) tmp(2) tmp(2)/tmp(1) tmp(3) tmp(4) tmp(4)/tmp(3) tmp(5) tmp(6) tmp(6)/tmp(5)]));
end;
tmp = quantile(triQ_,[0.25 0.75 0.01 0.99 0 1]); disp(sprintf('vol.s: 50 %% in [%0.1f %0.1f]%0.1f, 98 %% in [%0.1f %0.1f]%0.1f, all in [%0.1f %0.1f]%0.1f',[tmp(1) tmp(2) tmp(2)/tmp(1) tmp(3) tmp(4) tmp(4)/tmp(3) tmp(5) tmp(6) tmp(6)/tmp(5)]));
