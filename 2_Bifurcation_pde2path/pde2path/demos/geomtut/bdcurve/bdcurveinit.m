function p=bdcurveinit(nx,par) % spherical cap, init 
p=stanparam(); p.sw.spcalc=0; p.sw.bifcheck=0; % set stanparam, overwrite some
p.fuha.sG=@sGbdcurve; p.fuha.outfu=@cmcbra; p.fuha.e2rs=@e2rsbdry; 
pde=diskpdeo2(1,nx,round(nx/2)); % disk preimage discretization, pde-object 
% not stored, only p.DBC, p.tri and p.X (generated below) used subsequently
p.np=pde.grid.nPoints; p.nu=p.np; p.nt=pde.grid.nElements;  % store dimensions
p.sol.xi=1/p.nu; p.nc.neq=1; p.sw.sfem=-1; p.bcsw=0; % simplest case, see sG
po=pde.grid.p; x=po(1,:)'; y=po(2,:)'; u=0*ones(p.np,1); % set ICs 
p.u=[u; par]; p.X=[x,y,0*x]; p.tri=pde.grid.t(1:3,:)'; % store initial X and tri 
p.DIR=pde.grid.e(1:2,:)'; p.idx=unique(p.DIR(:)); % edges, and points for BCs 
p.sw.Xcont=1; p.sw.jac=0; % flag the X continuation, here with numjac 
p=oosetfemops(p); p.plot.auxdict={'H','V','A'}; p.plot.pstyle=-1; % userplot on 
p.u(p.nu+2)=getV(p,p.u); p.u(p.nu+3)=getA(p,p.u); % get initial vol & area
p.sw.bifcheck=0; p.sw.verb=2; p.plot.bpcmp=8; p.file.smod=1; 
p.usrlam=[]; p.nc.lammin=-0.1; p.nc.lammax=5; p.sw.spcalc=1; 
dsm=0.05; p.sol.ds=dsm; p.nc.dsmax=dsm; p.sw.rlong=1; p.sw.nobdref=0; 
p.nc.ilam=4; p.nc.lammax=1.5;  p.nc.lammin=-1.5; p.fuha.e2rs=@e2rsbdry; 