function un=p3interpol(xn,yn,zn,u,x,y,z,p)
% p3interpol: used for 3D plots
%
%  un=p3interpol(xn,yn,zn,u,x,y,z)
%
% un=interp of u from [x,y,z] to [xn,yn,zn]
% interpolate u def. on x,y,z to new mesh given in xn,yn,zn 
xv=reshape(x, size(x,1)*size(x,2), 1);
yv=reshape(y, size(y,1)*size(y,2), 1);
zv=reshape(z, size(z,1)*size(z,2), 1);
uv=reshape(u, size(u,1)*size(u,2), 1);
un= griddata3 (xv, yv, zv, uv, xn, yn, zn); 
end
