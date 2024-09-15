function p=acinit(p,r,hmax,par) % ac2D 
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
pde=tripdeo(r,hmax); % h as argument
p.cl=0; % (convenience) switch for later comparison with classical FEM 
p.plot.auxdict={'c','lambda','gamma','d'}; p.plot.pstyle=2; p.plot.cm='cool'; 
p.pdeo=pde; % domain and mesh
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e;
p.u=zeros(p.np,1); p.u=[p.u; par']; p.nc.ilam=2; p=setfemops(p);   
p.nc.nsteps=1000; p.sw.foldcheck=1; p.sol.ds=0.1; p.nc.dsmax=0.2; p.nc.lammax=1; 
p.plot.axis='image'; 
x0=[0, sqrt(3)/2]; x1=[-0.5, 0]; x2=[0.5, 0]; % find indizes of 3 points near corners 
% (in fact, exact corners, and indizes 1 3 5 known here, but this can be generalized) 
po=getpte(p)'; 
[k,dist]=dsearchn(po,[x0;x1;x2]); plotsol(p,1,1,1); view(0,90); hold on; 
plot(po(k,1),po(k,2),'*r'); p.pin=k; % save corner indizes for branch output
