function q=newmesh(p)
% NEWMESH: inquires to build new mesh and interpolate (u,tau) to it 
%
%  q=newmesh(p)
% (q is pde2path2 problem structure)
%
% See also stanmesh, newbmesh
q=p; 
msel=0; msel=asknu('Delaunay (0) or point-mesh(1)?',msel);
if(msel==0)
    xmin=min(p.mesh.p(1,:)); xmax=max(p.mesh.p(2,:)); 
    ymin=min(p.mesh.p(2,:)); ymax=max(p.mesh.p(2,:)); 
    hmax=min((xmax-xmin)/20,(ymax-ymin)/20);
    hmax=asknu('hmax:', hmax); q=stanmesh(q,hmax);
else nx=20;ny=20; nx=asknu('nx', nx);ny=asknu('ny', ny);
    q=stanmesh(q,nx,ny);
end
x=p.mesh.p(1,:); y=p.mesh.p(2,:); xn=q.mesh.p(1,:);yn=q.mesh.p(2,:);
u=p.mat.fill*p.u(1:p.nu);
for i=1:p.nc.neq % interpolate u to new mesh 
    z = u((i-1)*p.np+1:i*p.np); 
    ui=p2interpol(xn,yn,z,x,y); un((i-1)*q.np+1:i*q.np)=ui; 
end 
q.u=[un'; p.u(p.nu+1:end)]; q.nu = q.np*q.nc.neq;
if(q.sw.bcper~=0), q=rec2per(q,q.sw.bcper); end
plotsol(q,q.plot.ifig,q.plot.pcmp,q.plot.pstyle);  xi=1/q.np; 
if(q.sw.isw>1); fprintf('old xi=%g, new xi=%g\n',q.xi,xi); xi=asknu('xi:',q.xi); end 
q.xi=xi; 
otau=p.mat.fill*p.tau(1:p.nu); tau=zeros(q.np*q.nc.neq,1); 
for i=1:q.nc.neq % interpolate tau
  z = otau((i-1)*p.np+1:i*p.np);
  taui=p2interpol(xn,yn,z,x,y); 
  tau((i-1)*q.np+1:i*q.np) = taui; 
end
tau=[q.mat.drop*tau; p.tau(p.nu+1:end)]; % append aux variables to tau 
q.tau=tau/xinorm(tau,q.sol.xi,q.nc.nq,q.sol.xiq);


