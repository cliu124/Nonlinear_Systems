function uxy = tri2grid(p,t,u,x,y)
% tri2grip: interpolate (2D) u from triangulation to grid (obsolete ?)
xo=p(1,:); yo=p(2,:);  [X,Y]=meshgrid(x,y); 
[xi,yi,uxy]=griddata(xo,yo,u,X,Y); 
