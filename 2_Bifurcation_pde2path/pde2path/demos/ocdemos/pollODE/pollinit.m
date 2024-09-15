function p=pollinit(p,par) % init-routine for pollution demo 
p=stanparam(p); p.nc.neq=4; p.sw.jac=0; % numerical Jac 
p.fuha.sG=@pollsG; p.fuha.jcf=@polljcf; % rhs, objective value, 
p.fuha.outfu=@pollbra; % customized output (including objective function(s)) 
p.np=1; p.nu=4; p.sw.sfem=-1; p.vol=1; p.nc.neig=4; p=oosetfemops(p); 
p.sw.spcalc=0;  p.file.smod=10; p.sol.xi=0.005/p.np; p.plot.pstyle=5;
p.nc.dsmin=1e-6; p.nc.lammin=1e-10; p.nc.nsteps=150; p.plot.bpcmp=1; 
l1=-1; l2=-(par(2)+par(1)); % homogen. branch as initial guess 
z2=0.5*(1+par(1)-par(3)/(par(1)+par(2))); 
z1=z2*(1-z2); ov=ones(p.np,1);  u0=[z1*ov; z2*ov; l1*ov; l2*ov]; 
p.u=[u0; par']; p.nc.ilam=1; p.nc.imax=10; 
r=resi(p,p.u); fprintf('inires=%g\n',norm(r,'inf')); 
[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res); 
plotsol(p,1,1,p.plot.pstyle); 