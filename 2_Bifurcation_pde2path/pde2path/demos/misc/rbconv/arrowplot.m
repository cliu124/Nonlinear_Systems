function arrowplot(p,wnr,cnr)
% overlay density plot with arrow plot for vector field 
% from streamfunction in p
figure(wnr);n0=(cnr-1)*p.np+1; n1=cnr*p.np; 
u=p.mat.fill*p.u(1:p.nu);
pdeplot(p.mesh.p,p.mesh.e,p.mesh.t,'xydata',u(n0:n1));

% plot vector field arrows into figure
nx=20; ny=10; x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; 
xmin=min(x); xmax=max(x); ymin=min(y); ymax=max(y);
xg=linspace(xmin,xmax,nx);yg=linspace(ymin,ymax,ny);
ug=tri2grid(p.mesh.p,p.mesh.t,u(1:p.np),xg,yg);
[DX,DY] = gradient(ug);

figure(wnr); hold on;
quiver(xg,yg,-DY,DX,1,'LineWidth',1,'Color','black'); 
axis tight; hold off;