function p=sphinit(p,par,sw) % init for spherical vesicle 
p=stanparam(p); p=stanparamX(p); p.sw.Xcont=1; p.sw.jac=0;  % numJac
p.nc.neq=2;  p.file.smod=1; p.sw.verb=2; 
p.nc.neig=30; p.nc.eigref=-50; p.sw.bifcheck=2; p.sw.foldcheck=1; 
al=par(1); % used as radius 
[X,t]=subdivided_sphere(sw,'Radius',al); p.DBC=[]; p.idx=[]; % no BCs 
p.tri=[t(:,2) t(:,1) t(:,3)];p.X=X; p.Xold=X; 
[p.np,~]=size(p.X); p.nu=p.nc.neq*p.np; [p.nt,~]=size(p.tri); p.sol.xi=1/p.nu; 
u=zeros(p.np,1); u2=u+1/al; % initial H  
p=oosetfemops(p); p.u=[u;u2;par]; p.up=p.u;
p.plot.auxdict={'al','l1','l2','c0','A0','V0'}; p.plot.bpcmp=15; 
p.A0=getA(p,p.u); p.u(p.nu+5)=p.A0; p.V0=getV(p,p.u); p.u(p.nu+6)=p.V0;
p.fuha.outfu=@hvesbra; p.fuha.ufu=@refufu; p.nc.delbound=10; 
p.nc.sigr=0.02; p.nc.sigc=0.01; p.nc.Ab=max(doublearea(p.X,p.tri)/2); 
p.fuha.e2rs=@e2rsA; p.sw.rlong=1; p.nc.dsmax=1; p.sw.ips=2; p.nc.areafac=1;
p.nc.mu1=5; p.nc.mu2=0.1; p.nc.almin=0.1; p.nc.bisecmax=12; p.sw.para=1; p.DIR=[]; 
p.nc.foldtol=0.05; p.sw.foldcheck=0; p.nc.ilam=3; p.sol.ds=1;