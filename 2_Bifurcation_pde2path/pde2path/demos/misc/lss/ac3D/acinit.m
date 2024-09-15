function p=acinit(p,lx,ly,lz,nx,par) 
p=stanparam(p); screenlayout(p); p.nc.neq=1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
%pde=stanpdeo3D(lx,ly,lz,2*lx/nx);  % domain and mesh
pde=stanpdeo3D(lx,ly,lz, nx, round(ly/lx*nx), round(lz/lx*nx)); % alternate syntax 
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
p.u=zeros(p.np,1);p.u=[p.u; par']; p=setfemops(p);   
p.nc.nsteps=20;  p.sw.foldcheck=1; p.plot.auxdict={'c','lambda','gamma','d'}; 
p.plot.pstyle=1;  p.plot.cm='cool'; % plotting stuff
p.plot.levc={'blue','red'}; p.plot.alpha=0.4; p.plot.ng=20;