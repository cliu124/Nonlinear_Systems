function p=capinit(p,sw,par) % Helfrich caps, init 
p=stanparam(p); p=stanparamX(p); p.sw.jac=0;  % numJac
p.nc.neq=2; p.nc.neig=15; p.file.smod=1; p.sw.verb=2; 
p.nc.lammax=10^3; p.nc.lammin=-10^3; 
p.sw.bifcheck=2; p.sw.foldcheck=1; p.sf=1e0; al=par(1);
[p.X,p.tri]=subdiv_hemsphere(sw,'Base','icosahedron','Radius',al); 
p.tri=orient_outward(p.X,p.tri); p.tri=[p.tri(:,2) p.tri(:,1) p.tri(:,3)];
[p.np,~]=size(p.X); p.nu=p.nc.neq*p.np; [p.nt,~]=size(p.tri); u=zeros(p.np,1); 
p.sol.xi=1/p.nu; 
[q,~,~]=boundary_faces(p.tri); p.DIR=q; p.idx=unique([q(:,1);q(:,2)]); 
p=oosetfemops(p); u2=0*u+1/par(1); p.u=[u;u2;par]; p.up=p.u; p=retrigX(p); % 
p.plot.auxdict={'al','l1','c0','b'}; p.plot.bpcmp=8; 
p.nc.nq=0; p.fuha.qf=@qfV; p.fuha.qfder=@qjacV; p.sw.qjac=0; % for (later) PC 
p.fuha.outfu=@hcapbra; p.fuha.e2rs=@e2rs; 
p.nc.mu1=5; p.nc.mu2=0.2; p.nc.bisecmax=5; 
p.nc.foldtol=0.01; p.sw.foldcheck=0; p.nc.ilam=4; p.sol.ds=0.1; p.sw.ips=2; 