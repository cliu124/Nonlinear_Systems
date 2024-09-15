function un=p2interpol(xn,yn,u,x,y,p) 
% P2INTERPOL: interpolate u to new mesh.
%
%  un=p2interpol(xn,yn,u,x,y)
%
% * u,x,y: defined on mesh x,y  
% * xn,yn: new mesh 
xv=reshape(x, size(x,1)*size(x,2), 1);
yv=reshape(y, size(y,1)*size(y,2), 1);
uv=reshape(u, size(u,1)*size(u,2), 1);
un= griddata (xv, yv, uv, xn, yn); 
end
