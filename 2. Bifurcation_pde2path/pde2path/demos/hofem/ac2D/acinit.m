function p=acinit(p,lx,ly,nx,par,hofem) % ac2D 
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; p.hofem=hofem; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; p.nc.sf=1e3; 
pde=stanpdeo2D(lx,ly,2*lx/nx); % 2D rectangle pde object, h as argument
p.pdeo=pde; p=tri2six(p); % convert 3-node triang.to 6-node-triang: 
% add edge-midpoints, store 6-node-elems in p.hofem.tri, mesh in p.pdeo.grid 
p.pdeo.grid=setidssq(p.pdeo.grid); % set ids for boundary segments
%for i=1:4; identifyBoundarySegment(p.pdeo.grid,i); pause; end 
p.plot.auxdict={'c','lambda','gamma','d'}; p.plot.pstyle=1; p.plot.cm='cool'; 
%p.pdeo=pde; % domain and mesh
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e;
p.u=zeros(p.np,1); p.u=[p.u; par']; p.nc.ilam=2; p=setfemops(p);   
p.nc.nsteps=10; p.sw.foldcheck=1; p.sol.ds=0.1; p.nc.dsmax=0.2; p.nc.lammax=1; 
p.plot.axis='tight'; 
