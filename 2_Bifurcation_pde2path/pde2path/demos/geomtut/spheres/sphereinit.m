function p=sphereinit(par,r0,sw) % sphere 
p=stanparam(); p=stanparamX(p); p.sw.spcalc=0; p.sw.bifcheck=0; % set stanparam, overwrite some
[X,t]=subdivided_sphere(sw); X=r0*X; p.DBC=[]; p.idx=[];
p.tri=t; p.np=size(X,1); p.nu=p.np; p.nt=size(t,1);  % store dimensions
p.sol.xi=1/p.nu; p.nc.neq=1; u=0*ones(p.np,1); % set ICs 
p.u=[u; par]; p.X=X; p.tau=p.u; 
p=oosetfemops(p); p.plot.auxdict={'H','V','A'}; 
p.u(p.nu+2)=getV(p,p.u); p.u(p.nu+3)=getA(p,p.u); %  initial Vol and area 
p.u(p.nu+1)=-1/r0; p.sol.ds=0.1; p.nc.dsmax=4; p.nc.dlammax=4;
p.nc.ilam=[2 1]; p.nc.nq=1; p.fuha.qf=@qfV; p.fuha.qfder=@qjacV; 
p.plot.bpcmp=1; p.file.smod=2; p.sw.jac=0; 