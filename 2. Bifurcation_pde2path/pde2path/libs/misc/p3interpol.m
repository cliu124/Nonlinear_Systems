function un=p3interpol(xn,yn,zn,u,x,y,z,p)
% p3interpol: used for 3D plots
%
%  un=p3interpol(xn,yn,zn,u,x,y,z)
%
% un=interp of u from [x,y,z] to [xn,yn,zn]
% interpolate u def. on x,y,z to new mesh given in xn,yn,zn 
% behavior depends on p.sw.ips; 
% 0 (default): linear (interp)-linear (extrapol)
% 1: linear-nearest,    2: nearest-nearest 
ips=0; 
if nargin==8; try; ips=p.sw.ips; catch; end; end 
xv=reshape(x, size(x,1)*size(x,2), 1);
yv=reshape(y, size(y,1)*size(y,2), 1);
zv=reshape(z, size(z,1)*size(z,2), 1);
uv=reshape(u, size(u,1)*size(u,2), 1);
switch ips; 
  case 0; F=scatteredInterpolant(xv,yv,zv,uv,'linear','linear');
  case 3; F=scatteredInterpolant(xv,yv,zv,uv,'linear','nearest');
  case 2; F=scatteredInterpolant(xv,yv,zv,uv,'nearest','nearest');
end
un=F(xn,yn,zn);
if any(isnan(un)); fprintf('NaN in p2interpol, using nearest\n'); 
   F=scatteredInterpolant(xv,yv,zv,uv,'nearest','nearest');
   un=F(xn,yn,zn);
end