function p=animalinit(p) 
% init-routine 
p=stanparam(p); p.pstyle=2; 
p.nc.neq=2; p.fuha.G=@chemG; p.fuha.outfu=@chembra; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% geometry data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p.mesh.geo,p.mesh.bc]=animalgeo(0.025); hmax=0.2; %hmax=0.15; 
p=stanmesh(p,hmax); p.fuha.bc=@(p,u) p.mesh.bc; 
p.fuha.bcjac=@(p,u) p.mesh.bc;
p=setbmesh(p);
p.nc.dsmin=0.00001; p.nc.dlammax=0.25; p.nc.nsteps=50; p.nc.neig=min(20,p.np); 
p.sw.jac=0; p.sw.newt=1; p.nc.bisecmax=5; p.nc.amod=0; p.sw.bifcheck=1; p.plot.pcomp=1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% starting point %%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.nc.dlammax=0.05; p.sol.ds=0.05; u=1*ones(p.np,1); v=0.5*ones(p.np,1); 
p.nc.lammax=13; p.u=[u; v; 11.8];
p.nc.ilam=1;
p.tau=zeros(p.nc.neq*p.np+1,1); p.tau(1)=1; %plotsol(p,p.pfig,1,p.pstyle); 
p.vol=triint(ones(p.np,1),p.mesh.p,p.mesh.t); % needed in chembra 
p.sol.xi=1/p.np; pre=sprintf('%s',inputname(1)); p=setfn(p,pre); p.plot.cm='copper';
