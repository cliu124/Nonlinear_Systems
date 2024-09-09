function q=pa2pb(p,q) 
% pa2pb: interpolate u from p.pdeo.grid to q.pdeo.grid 
% 
%  q=pa2pb(p,q) 
% 
% poor man's adaption, intended just for comparison on two fixed grids
par=p.u(p.nu+1:end); u=p.u(1:p.nu); tau=p.tau(1:p.nu); 
uf=p.mat.fill*u; tauf=p.mat.fill*tau; lamd=p.tau(p.nu+1:end);   
onp=p.np; np=q.np;  ui=zeros(p.nc.neq*np,1); taui=ui; oldp=getpte(p); newp=getpte(q); 
switch size(oldp,1);
    case 2; xo=oldp(1,:); yo=oldp(2,:); xn=newp(1,:); yn=newp(2,:); 
for i=1:p.nc.neq
  ui((i-1)*np+1:i*np)=p2interpol(xn,yn,uf((i-1)*onp+1:i*onp),xo,yo); 
  taui((i-1)*np+1:i*np)=p2interpol(xn,yn,tauf((i-1)*onp+1:i*onp),xo,yo);
end
    case 3;
xo=oldp(1,:); yo=oldp(2,:); zo=oldp(3,:); xn=newp(1,:); yn=newp(2,:); zn=newp(3,:); 
for i=1:p.nc.neq
  ui((i-1)*np+1:i*np)=p3interpol(xn,yn,zn,uf((i-1)*onp+1:i*onp),xo,yo,zo); 
  taui((i-1)*np+1:i*np)=p3interpol(xn,yn,zn,tauf((i-1)*onp+1:i*onp),xo,yo,zo); 
end
end
q.u=[ui;par]; q.tau=[q.mat.drop*taui;lamd]; 
plotsol(p); figure(6); clf; plotsol(q,6,1,p.plot.pstyle); axis(q.plot.axis); 
r=resi(q,q.u); fprintf('inires=%g\n',norm(r,Inf));   imaxs=q.nc.imax; q.nc.imax=2*imaxs; 
[un,r,iter,Gu,Glam,q]=nloop(q,q.u); fprintf('res=%g\n',norm(r,q.sw.norm));
ud=un(1:q.nu)-ui(1:q.nu); maxd=max(abs(ud)); maxu=max(abs(un));
fprintf('max|u-ui|/max|u|=%g\n',maxd/maxu); q.nc.imax=imaxs; 
