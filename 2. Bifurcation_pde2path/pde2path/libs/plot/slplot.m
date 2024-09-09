function slplot(p,ng,u,fs) 
%  slplot: slice plot (3D) 
%
%  slplot(o,ng,u,fs)  
%
% o=pde-object, ng=#grid points for interp., u=sol, fs=fontsize
%  rem call fig-nr before sl-plot cause directly used in hoplot
gp=getpte(p); 
x1=min(gp(1,:)); x2=max(gp(1,:)); y1=min(gp(2,:)); y2=max(gp(2,:)); 
z1=min(gp(3,:)); z2=max(gp(3,:)); 
xv=linspace(x1,x2,ng); yv=linspace(y1,y2,ng); zv=linspace(z1,z2,ng); 
[X,Y,Z]=meshgrid(xv,yv,zv);
up=p3interpol(X,Y,Z,u,gp(1,:),gp(2,:),gp(3,:),p);
xsl=(x1+x2)/2; ysl=(y1+y2)/2; zsl=[(z1+z2)/2,0];
%figure(fnr); clf; 
slice(X,Y,Z,up,xsl,ysl,zsl);
view(3); axis([x1 x2 y1 y2 z1 z2]); 
xlabel('x','fontsize',fs); ylabel('y','fontsize',fs); 
zlabel('z','fontsize',fs); box on;set(gca,'FontSize',fs); 
