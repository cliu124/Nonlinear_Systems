function [xy,tri,bndmesh,options,geomfunc] = mesh_rect3D(N,options,ustruct)
if ustruct ~= 0
    x = rand(1,N^3); y = rand(1,N^3); z = rand(1,N^3);
    %add edges
    x = [x 0 0 0 0 1 1 1 1];
    y = [y 0 0 1 1 0 0 1 1];
    z = [z 0 1 0 1 0 1 0 1];
    x = [x rand(1,(N-1)*4)];
    y = [y ones(1,N-1) ones(1,N-1)  zeros(1,N-1) zeros(1,N-1)];
    z = [z ones(1,N-1) zeros(1,N-1) zeros(1,N-1) ones(1,N-1)];
    x = [x ones(1,N-1) ones(1,N-1)  zeros(1,N-1) zeros(1,N-1)];
    y = [y rand(1,(N-1)*4)];
    z = [z ones(1,N-1) zeros(1,N-1) zeros(1,N-1) ones(1,N-1)];
    x = [x ones(1,N-1) ones(1,N-1)  zeros(1,N-1) zeros(1,N-1)]';
    y = [y ones(1,N-1) zeros(1,N-1) zeros(1,N-1) ones(1,N-1)]';
    z = [z rand(1,(N-1)*4)]';
    if ustruct == 2
      %%add faces
      %[XYZ,YZX] = meshgrid(1/N:1/N:1-1/N); XYZ = repmat(XYZ(:),2,1); YZX = repmat(YZX(:),2,1); XYZ0 = [zeros((N-1)^2,1); ones((N-1)^2,1)];
      XYZ = rand((N-1)^2,1); YZX = rand((N-1)^2,1); XYZ = repmat(XYZ(:),2,1); YZX = repmat(YZX(:),2,1); XYZ0 = [zeros((N-1)^2,1); ones((N-1)^2,1)];
      x = [x; XYZ0; YZX; XYZ]; 
      y = [y; XYZ; XYZ0; YZX]; 
      z = [z; YZX; XYZ; XYZ0];
    end;
    %tri = delaunay3(x,y,z);
    tri=delaunayn([x y z]); 
    tri = elem_fix(tri,[x y z]);
else
[x,y,z] = meshgrid(0:1/N:1); x=x(:); y=y(:); z=z(:);
[xn,yn,zn] = meshgrid(1:N); xn=xn(:); yn=yn(:); zn=zn(:);
trin = yn+(N+1)*(xn-1)+(N+1)^2*(zn-1);
trin = [trin trin+1 trin+N+1 trin+N+2];
trin = [trin trin+(N+1)^2];
tri = [trin(:,[1 5 7 6]); trin(:,[8 2 4 3]); trin(:,[2 7 6 1]); trin(:,[2 7 1 3]); trin(:,[2 7 8 6]); trin(:,[2 7 3 8])];
end;
[C,I] = sort(rand(size(tri,1),1)); tri = tri(I,:); %shuffle elements
[C,I] = sort(rand(size(x))); x = x(I); y = y(I); z = z(I); %shuffle nodes
[C,I2] = sort(I); tri = I2(tri);

xy = [x y z];
bndmesh = [];
options.area = 1;
if nargout == 5
	geomfunc = [];
end;