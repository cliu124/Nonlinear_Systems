function p=acgcinit(p,par,lx,nx,dim,brute) % init-routine for acgc 
p=stanparam(p); p.nc.neq=1; p.fuha.outfu=@acgcbra; p.jfac=0; p.sw.verb=2; p.nc.mu1=1; 
if brute;  p.jfac=1; % brute force (full Jac) for global terms 
else p.fuha.lss=@gclss; p.fuha.blss=@gcblss; % Sherman-Morrison version 
     p.sw.eigssol=1; % use SM also in eigs
end
p.sw.runpar=0;  % switch off parfor in pmnewtonloop (clashes with global vars) 
switch dim 
    case 1;  p.nc.neig=10; p.plot.pstyle=1; 
       pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; p.Om=2*lx; % domain and mesh
    case 2; ly=1.1*pi/2; ny=nx*round(ly/lx); p.nc.neig=30; p.plot.pstyle=3;
       pde=stanpdeo2D(lx,ly,nx,ny); p.pdeo=pde; p.Om=4*lx*ly; % domain and mesh
end
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t;  
p.mesh.be=e; p.mesh.nt=size(t,2); % background-mesh (for oomeshadac)
% PDE simple implementation
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
p.tau=1; p.sol.xi=1/p.np;  p.file.smod=5; 
p.sw.bifcheck=2; p.nc.lammin=-1; p.nc.lammax=6; p.plot.bpcmp=0;
p.nc.imax=10; p.nc.dsinciter=p.nc.imax/4; p.nc.dsmax=0.1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% starting point %%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.nc.ilam=2; po=getpte(p); x=po(1,:)'; p.u=0*x; p.u=[p.u;par]; 
p.sol.ds=0.01; p=setfemops(p); 
res=norm(resi(p,p.u), p.sw.norm); fprintf('initial res=%g\n',res);  
[p.u,res,iter,~,~]=nloop(p,p.u); plotsol(p,1,1,p.plot.pstyle);
fprintf('first res=%g with iter=%i\n',res,iter);