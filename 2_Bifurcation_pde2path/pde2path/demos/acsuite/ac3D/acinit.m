function p=acinit(p,lx,ly,lz,nx,par,varargin) % ac3D 
if nargin>6; sw=varargin{1}; else sw.sym=0; sw.ref=0; end 
p=stanparam(p); screenlayout(p); p.nc.neq=1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
pde=stanpdeo3D(lx,ly,lz,2*lx/nx,sw);  % domain and mesh
%pde=stanpdeo3D(lx,ly,lz,nx,round(ly/lx*nx),round(lz/lx*nx),sw); % alt.syntax 
p.plot.pstyle=1;  p.plot.cm='cool'; % plotting stuff
p.plot.levc={'blue','red'}; p.plot.alpha=0.1; p.plot.ng=20;
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
p.u=zeros(p.np,1);p.u=[p.u; par']; p=setfemops(p);   
p.nc.nsteps=20;  p.sw.foldcheck=1; p.plot.auxdict={'c','lambda','gamma','d'}; 