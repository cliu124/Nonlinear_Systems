function p=spcapinit(nx,par) % spherical cap, init, legacy setup 
p=stanparam(); p.sw.spcalc=0; p.sw.bifcheck=0; % set stanparam, overwrite some
pde=diskpdeo2(1,nx,round(nx/2)); % disk preimage discretization 
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.nt=pde.grid.nElements; 
p.sol.xi=1/p.nu; p.nc.neq=1; p.sw.sfem=-1; % store dimensions
p.fuha.outfu=@cmcbra; p.fuha.e2rs=@e2rs; % branch data, refinement selector
p.sw.Xcont=1; p.plot.pstyle=-1; % call userplot for plotting
po=getpte(p); x=po(1,:)'; y=po(2,:)'; u=0*ones(p.np,1); % initial soln 
p.u=[u; par]; p.X=[x,y,0*x]; % set IC, including X 
p=oosetfemops(p); % here constant mass and stifness matrices (until mesh changes) 
p.plot.auxdict={'H','V','alpha','A'}; p.u(p.nu+4)=getA(p,p.u); % initial area 
