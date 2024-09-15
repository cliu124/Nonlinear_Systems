clc; clear all; close all;
times = []; rand('state',0);
Nlist = 20; do3D = false; %40:15:100;  %%110:20:190;
%Nlist = 15; do3D = true;
for N=Nlist
    disp(sprintf('N=%0.0f',N));
    options = gen_options();
    options.qualM = 1;
    %options.qualP = -pi/2*1.2;
    options.prag_adapt = 2;
    options.debug = 0;
    options.RMedg = 0;
    options.qualRM = 0;
    options.fastRM = 2;
    options.consRM = 1;
    options.greedyCLR = 2;
    options.advRFN = 0;
    options.smpRFN = 0;
    options.verbose = 1;
    options.minA = eps;
    options.consRFN = 0;
    options.nosparse = true;
    options.innerit = 15;
    %options.prag_crs = true;
    options.mntn_bks = 1;
    %options.fastRFN = 0;
    rand('state',0);

    if do3D
    	options.consRM = 0;
    	options.smpRFN = 1;
    	options.consRFN = 2;
    	options.prag_crs = true;
    	options.swap3D = false;
    	options.innerit = 20;
    [xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(2*N,options,0);
    metric = @(x_) metric_sphere(x_, 0.2, 0.5, [eps eps 1-eps], 1/N, 1/N/10);	
    else
    [xy,tri,bndmesh,options,geomfunc] = mesh_rect(2*N,options);
    metric = @(x_) metric_circ(x_, 0.2, 0.5, [eps eps], 1/N, 1/N/10);
    end;
    Nmetric = metric(xy);
    %metric = @(x_) metric_shock(x_, 0.1, 0.1, repmat(0.5,1,size(xy,2)), 5*1/N, 1/N/10);
    %metric = @(x_) metric_uniform(x_, repmat(1/N,1,size(xy,2)));
    [C,I] = min(abs(xy(:,1)-0.25)); thex = xy(I,1);
    if size(xy,2) == 2 %2D set line x==0.25
    bndmesh = bndmesh_polygon(tri,xy,bndmesh,options);
    bndmesh.triID = ones(size(tri,1),1); trixy = (xy(tri(:,1),:)+xy(tri(:,2),:)+xy(tri(:,3),:))/size(tri,2); 
    edgbxy = (xy(bndmesh.edg(:,1),:)+xy(bndmesh.edg(:,2),:))/2.;
    [edg,edg2tri,tri2edg] = bks_all(tri); edgxy = (xy(edg(:,1),:)+xy(edg(:,2),:))/2.;
    for jj=1:2
    bndmesh.triID(trixy(:,jj) < thex) = bndmesh.triID(trixy(:,jj) < thex)+jj;
    Ileft = and(or(edgbxy(:,3-jj) < eps,1.-eps < edgbxy(:,3-jj)),edgbxy(:,jj)<thex);
    bndmesh.IDs(Ileft) = 5+jj;
    Icntr = and(thex-eps<edgxy(:,jj),edgxy(:,jj)<thex+eps);
    bndmesh.edg = [bndmesh.edg; edg(Icntr,:)];
    bndmesh.IDs = [bndmesh.IDs; repmat(13+jj,nnz(Icntr),1)];
    end;
    else %3D set face at x==0.25
    bndmesh = bndmesh_polyhedron(tri,xy,bndmesh,options);
    [fac,fac2tri,tri2fac] = bks_all3D(tri);
    bndmesh.triID = ones(size(tri,1),1); trixy = (xy(tri(:,1),:)+xy(tri(:,2),:)+xy(tri(:,3),:)+xy(tri(:,4),:))/size(tri,2); 
    facxyb = (xy(bndmesh.fac(:,1),:)+xy(bndmesh.fac(:,2),:)+xy(bndmesh.fac(:,3),:))/3.;
    facxy = (xy(fac(:,1),:)+xy(fac(:,2),:)+xy(fac(:,3),:))/3.;
    for jj=1:2
    bndmesh.triID(trixy(:,jj) < thex) = bndmesh.triID(trixy(:,jj) < thex)+round(0.5*jj*(jj-1))+1;
    Ileft = and(eps < facxyb(:,jj), facxyb(:,jj)<thex);
    bndmesh.IDs(Ileft) = bndmesh.IDs(Ileft)+5+10*jj;
    Icntr = and(thex-eps<facxy(:,jj),facxy(:,jj)<thex+eps);
    bndmesh.fac = [bndmesh.fac; fac(Icntr,:)];
    bndmesh.IDs = [bndmesh.IDs; repmat(13+(10*(jj-1))^2,nnz(Icntr),1)];
    if jj==2
    	facxyb = (xy(bndmesh.fac(:,1),:)+xy(bndmesh.fac(:,2),:)+xy(bndmesh.fac(:,3),:))/3.;
    	Ileft = or(and(and(thex-eps < facxyb(:,1), facxyb(:,1) < thex+eps), ...
    	thex-eps < facxyb(:,2)),and(thex-eps < facxyb(:,1), ...
    	and(thex-eps < facxyb(:,2), facxyb(:,2) < thex+eps)));
    	bndmesh.IDs(Ileft) = bndmesh.IDs(Ileft) + 1;
    end;
    end;
    numel(unique(bndmesh.IDs)) %11,20
    Nmetric = Nmetric*2.^(1./3.);
    end;
    tic;
    for i=1:options.outerit
	[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
    end;
    times(end+1) = toc/size(xy,1)*1e3
end;
if numel(Nlist) > 1
     plot(Nlist/Nlist(end)*size(xy,1)/1e3,times,'.b-'); ylabel('seconds/1k nodes'); xlabel('k nodes'); title(sprintf('pragmatic takes %0.2f seconds/1k nodes (for 250k)',24/250)); %1e3/times(end)/(250e3/24) 20% of pragmatic speed
else
    colors = {'r','b','g','k'};
    if size(xy,2) == 2
	trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','.'); xlim([-0.1 1.1]); ylim([-0.1 1.1]);
	trixy = [mean(reshape(xy(tri(:),1),size(tri)),2) mean(reshape(xy(tri(:),2),size(tri)),2)];
	hold on; for i=min(bndmesh.triID):max(bndmesh.triID)
		I = bndmesh.triID==i;
		plot(trixy(I,1),trixy(I,2),['s' colors{i+1-min(bndmesh.triID)}],'markersize',3);
	end; hold off;
    else
	trimesh(bndmesh.fac(bndmesh.IDs==13,:),xy(:,2),xy(:,3))
    end;
end;
trix = mean(reshape(xy(tri,1),size(tri)),2); [all(trix(bndmesh.triID==1) > thex) all(trix(bndmesh.triID==2) < thex)]