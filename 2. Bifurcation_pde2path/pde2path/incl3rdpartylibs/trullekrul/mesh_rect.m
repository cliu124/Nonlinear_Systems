function [xy,tri,bndmesh,options,geomfunc] = mesh_rect(N,options,ustruct)
if nargin ~= 2
x = rand(1,N*N);
y = rand(1,N*N);
x = [x 0 0 1 1 zeros(1,N) rand(1,N) ones(1,N) rand(1,N)]';
y = [y 0 1 1 0 rand(1,N) ones(1,N) rand(1,N) zeros(1,N)]';
tri = delaunay(x,y);
tri = elem_fix(tri,[x y]);
else
[x,y] = meshgrid(0:1/N:1); x=x(:); y=y(:);
[xn,yn] = meshgrid(1:N); xn=xn(:); yn=yn(:);
trin = yn+(N+1)*(xn-1); 
trin = [trin trin+1 trin+N+1 trin+N+2];
tri = [trin(:,[2 1 3]); trin(:,2:4)];
end;
[C,I] = sort(rand(size(tri,1),1)); tri = tri(I,:); %shuffle elements
[C,I] = sort(rand(size(x))); x = x(I); y = y(I); %shuffle nodes
[C,I2] = sort(I); tri = I2(tri);
xy = [x y];
bndmesh = [];
options.area = 1.;
if nargout == 5
	geomfunc = [];
	%geomfunc = @(x_) geom_rect(x_,1,1); 
end;