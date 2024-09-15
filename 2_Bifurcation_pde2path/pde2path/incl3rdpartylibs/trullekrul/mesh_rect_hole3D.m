function [xy,tri,bndmesh,options,geomfunc] = mesh_rect_hole3D(N,dxyz,r,options)
%N = 3; dxyz = [12 2 2]; r=1; options = gen_options();
[xy,tri,bndmesh_,options_,geomfunc_] = mesh_rect3D(3,options,0);
trixy = reshape(xy(tri(:),:),size(tri,1),12); trixy = [mean(trixy(:,1:4),2) mean(trixy(:,5:8),2) mean(trixy(:,9:12),2)];
I = all(abs(trixy(:,1:2)-0.5) < 1/6,2);
tri = tri(not(I),:); 
xy(:,1) = dxyz(1)*(xy(:,1)-0.5);
xy(:,2) = dxyz(2)*(xy(:,2)-0.5);
xy(:,3) = dxyz(3)*(xy(:,3)-0.5);
I = sqrt(sum(xy(:,1:2).^2,2)) < sqrt(dxyz(2)/3+(dxyz(1)/6)^2)+options.geomtol;
xy(I,1) = r*2^-0.5*sign(xy(I,1));
xy(I,2) = r*2^-0.5*sign(xy(I,2));

bndmesh = bndmesh_polyhedron(tri,xy,[],options);
facxy = reshape(xy(bndmesh.fac(:),:),size(bndmesh.fac,1),9); facxy = [mean(facxy(:,1:3),2) mean(facxy(:,4:6),2) mean(facxy(:,7:9),2)];

Nmetric = [ones(size(xy,1),1) zeros(size(xy,1),1) ones(size(xy,1),1) ...
	  zeros(size(xy,1),1) zeros(size(xy,1),1) ones(size(xy,1),1)]*N;

distcyl = @(x_) geom_cyl(x_,r);
distR   = @(x_) geom_rect(x_,dxyz,-dxyz/2);
geomfunc = @(x_) geom_diff(distR(x_),distcyl(x_));
options.area = 0;
options_ = options;
options_.debug = 2;
%options_.fstRFN = 1;
options_.fstRFN = 0;
options_.advRFN = 1;
options_.RMnd3D = 1;
options_.innerit = 20;
%adapt
[tri,xy,Nmetric,bndmesh] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options_);
%trimesh(tri',xy(:,1),xy(:,2),xy(:,3))
%shuffle nodes
[C,I] = sort(rand(size(xy,1),1));
xy = xy(I,:);
[C,I] = sort(I);
tri = I(tri);
bndmesh.fac = sort(I(bndmesh.fac),2);
