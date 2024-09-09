function [xy,tri,bndmesh,options,geomfunc] = mesh_sphere(N,options,ustruct)
N = round(N/2);
[xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(N,options,ustruct);
xy(:,1) = 1-xy(:,1);
%xy(:,2) = 1-xy(:,2);
xy(:,3) = 1-xy(:,3);
[tri,xy] = min_reflect(tri,xy,1);
[tri,xy] = min_reflect(tri,xy,2);
[tri,xy] = min_reflect(tri,xy,3);
%transform cube to sphere
bndmesh = bndmesh_polyhedron(tri,xy,[],options); x = xy(:,1); y = xy(:,2); z = xy(:,3);
%x = 2*(x - 0.5); y = 2*(y - 0.5); z = 2*(z - 0.5);
xy = [x.*(1.-0.5*(z.^2 + y.^2) + z.^2.*y.^2/3.).^(1/2) ...
      y.*(1.-0.5*(x.^2 + z.^2) + x.^2.*z.^2/3.).^(1/2) ...
      z.*(1.-0.5*(y.^2 + x.^2) + y.^2.*x.^2/3.).^(1/2)];
r = sqrt(sum(xy.^2,2));
phi = atan2(xy(:,2),xy(:,1));
theta = atan2(xy(:,3),sqrt(sum(xy(:,1:2).^2,2)));
r = sign(r).*abs(r).^2; %1.5;
xy = [r.*cos(phi).*cos(theta) r.*sin(phi).*cos(theta) r.*sin(theta)];
xy = xy/2 + 0.5;

bndmesh.crnds = []; bndmesh.IDs(:) = 1;
options.area = 0;
geomfunc = @(x_) geom_circ(x_,0.5,[0.5 0.5 0.5]);
geomfunc = {geomfunc,[]};

function [tri,xy] = min_reflect(tri,xy,dim)
good = xy(:,dim) ~= min(xy(:,dim));
newtri = tri(:,[2 1 3 4])+size(xy,1);
flip = ones(1,3); flip(dim) = -1;
xy = [xy; xy(good,:).*repmat(flip,sum(good),1)];
I = [true(size(good)); good];
old2new = zeros(size(I)); old2new(I) = 1:sum(I);
old2new(not(I)) = find(not(good));
tri = [tri; old2new(newtri)];