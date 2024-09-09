function p=acinit(p,lx,nx,par,dim,sw,del) % init with prep. of mesh-adaption 
p=stanparam(p); screenlayout(p); p.nc.neq=1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
switch dim; 
    case 1; pde=stanpdeo1Db(0,lx,lx/nx); p.plot.axis=[0 1 0 1]; 
    case 2; pde=diskpdeo(lx,nx,del); % refine all longest edges for more symm.      
            rlongs=pde.grid.rlong; pde.grid.rlong=1; idx=1:size(pde.grid.t,2);
            pde.grid.refineMesh(idx); pde.grid.rlong=rlongs;
    case 3; pde=stanpdeo2D(lx,lx/2,lx/nx); 
end
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu);
p.mesh.bp=pde.grid.p; p.mesh.bt=pde.grid.t; p.mesh.be=pde.grid.e; 
[po,t,e]=getpte(p);  x=po(1,:)'; x0=0.5; a=0.0001; b=16*(1-a); u0=a+b*(x-x0).^4; 
switch sw; 
    case 0; u0=0*u0; 
    case 1; u0=1+0*x; 
    case 2; y=po(2,:)'; r=x.^2+y.^2;  b=1-a; u0=a+b*r.^2; 
end 
p.u=[u0; par']; % initial guess, and parameters 
p.usrlam=[]; p=setfemops(p); p.mesh.maxt=100; 
p.plot.auxdict={'\lambda','\gamma'}; 
p.nc.ilam=1; p.nc.lammin=-1; p.nc.lammax=50; 
p.sol.ds=0.1; p.nc.dsmax=2; p.nc.sf=1000; p.nc.del=1e-2; 