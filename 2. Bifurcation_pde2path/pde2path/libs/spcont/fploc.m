function p=fploc(p) 
% fploc: localize foldpoint via extended system 
%
% p=fploc(p) 
q=p; u0=p.u; p.sol.ptype=2; p.u=[u0;0]; % add dummy parameter
p=spcontini(p,length(p.u)-p.nu,'',0); p.fuha.spjac=@spjac; 
r=resi(p,p.u);res=norm(r,p.sw.norm); fprintf('ini-res=%g\n',res); 
[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res);   
p.sw.spcont=0; p.nu=p.nu/2; % mimic spconexit here 
p.nc.neq=p.nc.neq/2; p.nc.nq=(p.nc.nq-1)/2; % set regular case sizes
p.u=[p.u(1:p.nu);p.u(2*p.nu+1:2*p.nu+p.naux)]; % non-fold/branch point variables
p.nc.ilam=p.nc.ilam(2:2+p.nc.nq); % remove current primary parameter
u=p.u(1:end-1); % take dummy parameter away and store p.u
p=q; p.u=u; % set all parameters to the original ones, set p.u to u
r=resi(p,u); [Gu,Glam]=getder(p,u,r); % compute new tangent
amat=genamat(p,Gu,Glam,p.tau,p.sol.xi,p.sol.xiq); 
[tau,p]=p.fuha.blss(amat,[zeros(p.nu+p.nc.nq,1);1],p); 
tau=tau/xinorm(tau,p.sol.xi,p.nc.nq,p.sol.xiq);  p.tau=tau; 
p.branch=[p.branch [bradat(p); p.fuha.outfu(p,p.u)]]; % add branch data 
p.sol.ptype=1; p.fuha.savefu(p); % save 
p.file.bcount=p.file.bcount+1; p.file.count=p.file.count+1; % adjust counters