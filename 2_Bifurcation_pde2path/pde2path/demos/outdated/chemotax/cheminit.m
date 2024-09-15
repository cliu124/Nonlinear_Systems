function p=cheminit(p) 
% init-routine 
p=stanparam(p); pre=sprintf('%s',inputname(1)); p=setfn(p,pre);
p.nc.neq=2; p.fuha.G=@chemG; p.fuha.outfu=@chembra; 
lx=0.5; ly=2; p.mesh.geo=rec(lx,ly); p=stanmesh(p,0.075); 
bc=gnbc(p.nc.neq,4,zeros(p.nc.neq),zeros(p.nc.neq,1)); 
p.fuha.bc=@(p,u) bc; p.fuha.bcjac=@(p,u) bc;
p.nc.dsmin=0.00001; p.nc.dlammax=0.5; p.nc.nsteps=20; p.nc.neig=min(50,p.np); 
p.sw.spcalc=0; p.sw.jac=0;  % use numjac! 
p.sw.newt=1; p.nc.bisecmax=15; p.nc.amod=0; p.sw.bifcheck=1; 
p.plot.pcmp=2; p.sw.para=2; p.sol.ds=0.1;
p.nc.dlammax=2; p.nc.lammax=25; u=1*ones(p.np,1); v=0.5*ones(p.np,1); 
p.u=[u; v; 10]; p.nc.ilam=1;
plotsol(p,p.plot.pfig,1,p.plot.pstyle);
p.vol=triint(ones(p.np,1),p.mesh.p,p.mesh.t); % needed in chembra 
p.sol.xi=1/p.np; p.mesh.maxt=2*p.mesh.nt; 



