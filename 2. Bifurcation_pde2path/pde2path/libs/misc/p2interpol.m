function un=p2interpol(xn,yn,u,x,y,p) 
% P2INTERPOL: interpolate u to new mesh.
%
%  un=p2interpol(xn,yn,u,x,y)
%
% * u,x,y: defined on mesh x,y  
% * xn,yn: new mesh 
% * matlab interface, see octcomp/p2interpol for octave version 
%
% behavior depends on p.sw.ips; 
% 0 (default): linear (interp)-linear (extrapol)
% 1: linear-nearest,    2: nearest-nearest 
ips=0; 
if nargin==6; try; ips=p.sw.ips; catch; end; end 
xv=reshape(x, size(x,1)*size(x,2), 1);
yv=reshape(y, size(y,1)*size(y,2), 1);
zv=reshape(u, size(u,1)*size(u,2), 1);
switch ips; 
  case 0; F=scatteredInterpolant(xv,yv,zv,'linear','linear');
  case 1; F=scatteredInterpolant(xv,yv,zv,'linear','nearest');
  case 2; F=scatteredInterpolant(xv,yv,zv,'nearest','nearest');
end
un=F(xn,yn);
if any(isnan(un)); fprintf('NaN in p2interpol, using nearest\n'); 
   F=scatteredInterpolant(xv,yv,zv,'nearest','nearest');
   un=F(xn,yn,zn);
end
