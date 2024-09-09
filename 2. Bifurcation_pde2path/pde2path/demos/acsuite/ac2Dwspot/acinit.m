function p=acinit(p,lx,ly,nx,par,varargin) %  
if nargin>5; sw=varargin{1}; else sw.sym=0; sw.ref=0; end 
p=stanparam(p); screenlayout(p); p.nc.neq=1; 
p.sw.sfem=-1; p.fuha.sG=@sGws; p.fuha.sGjac=@sGwsjac; 
pde=stanpdeo2D(lx,ly,nx,round(ly/lx*nx),sw); % alternate syntax 
p.plot.pstyle=1;  p.plot.cm='cool'; % plotting stuff
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
p.u=zeros(p.np,1);p.u=[p.u; par']; p=setfemops(p);   
p.nc.nsteps=20;  p.sw.foldcheck=1; p.plot.auxdict={'c','lambda','gamma','xi'}; 
p.plot.view=[5 40];  