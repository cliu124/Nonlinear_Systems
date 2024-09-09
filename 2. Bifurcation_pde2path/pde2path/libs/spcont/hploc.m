function p=hploc(p,varargin) 
% hploc: localize Hopf point via extended system
%
% p=hploc(p) 
% p=hploc(p,eigref)   (to choose im-part of desired eigenvalue) 
eigref=1; if nargin>1; eigref=varargin{1}; end 
p.sol.ptype=3; p.u=[p.u;0]; % add dummy parameter
aux.eigref=eigref; p=hpcontini(p,length(p.u)-p.nu,'',0,aux); 
r=resi(p,p.u); res=norm(r,p.sw.norm); fprintf('ini-res=%g\n',res); 
[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res);   
p.sw.spcont=0; p.nu=(p.nu-1)/3; % mimic hpcontexit here 
p.nc.neq=p.nc.neq/3; p.nc.nq=(p.nc.nq-1)/2; % set regular case sizes
p.u=[p.u(1:p.nu);p.u(3*p.nu+2:3*p.nu+p.naux)]; % non-HP point variables
p.nc.ilam=p.nc.ilam(2:2+p.nc.nq); % remove current primary parameter
p.tau=[p.tau(1:p.nu);p.tau(3*p.nu+2:3*p.nu+1+p.nc.nq);0];
p.branch=[p.branch [bradat(p); p.fuha.outfu(p,p.u)]]; % add branch data 
p.sol.j0=eigref; p.fuha.savefu(p); % save 
p.file.bcount=p.file.bcount+1; p.file.count=p.file.count+1; % adjust counters