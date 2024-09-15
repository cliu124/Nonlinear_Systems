function qual = elem_cond(area,v1,v2,v3)
gdim = nargin-1;
inds = reshape(1:gdim ^2,[gdim gdim]);
G = zeros(size(v1,1),gdim^2);
if gdim == 2
J = [v1 v2];
%area = area*(sqrt(3)/4)/0.5;
else
J = [v1 v2 v3];
%area = area*(sqrt(2)/12)/(1/6);
end;
for i1=1:gdim
for i2=1:gdim
for i3=1:gdim
G(:,inds(i1,i3)) = G(:,inds(i1,i3)) + J(:,inds(i1,i2)).*J(:,inds(i2,i3));
end;
end;
end;
detJ = area;

%2D
xy1 = rand(2,1);
xy2 = rand(2,1);
metric = rand(2,2); metric = metric*metric';
metric2 = metric*metric;
xy1_ = metric*xy1;
xy2_ = metric*xy2;
Lsq_ = sum(xy1_.^2)
A_ = det([xy1_ xy2_])
%using metric with unit of inverse squared length
Lsq = xy1'*metric2*xy1
A = det([xy1 xy2])*sqrt(det(metric2))
%3D
xy1 = rand(3,1);
xy2 = rand(3,1);
xy3 = rand(3,1);
metric = rand(3,3); metric = metric*metric';
metric2 = metric*metric;
xy1_ = metric*xy1;
xy2_ = metric*xy2;
xy3_ = metric*xy3;
Lsq_ = sum(xy1_.^2)
A_ = sum(cross(xy1_,xy2_).^2)
V_ = det([xy1_ xy2_ xy3_])
%using metric with unit of inverse squared length
Lsq = xy1'*metric2*xy1
A = det([xy1 xy2]'*metric2*[xy1 xy2])
V = det([xy1 xy2 xy3])*sqrt(det(metric2))


%stabilisation
A = rand(2,2);
A = A*A';
[v,L] = eig(A); 
L = diag([3. 1.]); A = v*L*v';
t = linspace(0,2*pi,100);
xy = repmat(v(:,1),1,numel(t)).*repmat(cos(t),2,1)*L(1,1)+repmat(v(:,2),1,numel(t)).*repmat(sin(t),2,1)*L(2,2);
plot(xy(1,:),xy(2,:),'-r',[0 v(1,1)*L(1,1)],[0 v(2,1)*L(1,1)],'-k',[0 v(2,1)*L(2,2)],[0 v(2,2)*L(2,2)],'-k'); axis equal;

AR = max(L(1,1)/L(2,2),L(2,2)/L(1,1))
n = rand(2,1);
Ai = inv(A);
ni = Ai*n;
niL = 1/sqrt(sum(ni.^2));
hold on; plot([0 n(1)*niL],[0 n(2)*niL],'-b'); hold off;
niL2 = 1/sqrt(n'*Ai*Ai*n)

%load for_debug3D.mat; geomfunc = []; Nmetric = metric_sqrt(Nmetric); tri = elem_fix(tri,xy);
% [tri,xy,Nmetric,bndmesh] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc ,options);