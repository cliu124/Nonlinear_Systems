%% demo schnaktravel
close all; keep pphome;
%% init and findbif 
p=[]; p=schnak_init(p,1,40); 
p.sol.ds=-0.01; p.nc.dsmin=0.001; p.nc.dsmax=0.1; 
p.nc.amod=0; p.nc.nsteps=10;  p=findbif(p,3); 
%% swibra to stripes 
p=swibra('p','bpt1','q',-0.1); p.nc.lammin=1.3; p.nc.nsteps=60; p=cont(p);
%% continue in domain size: par(3) multiplies diffusion matrix
clf(2); p=swiparf('q','pt30','q2',3); p.nc.lammin=0.1;
p.sol.ds=-1e-2; p.nc.dsmax=5e-2; p.file.smod=10;p=cont(p,45);
%% continue from symmetry breaking bif with Neumann bc
% (to contrast with pBC below) 
p=swibra('q2','bpt2','q3',0.1); p.nc.nsteps=10; p=cont(p);
%% switch to cylinder and add phase condition
% change to cylinder geometry, init with last point from NBC run in par(3)
p=loadp('q2','pt40','r',1); p=resetc(p); aux=p.u(p.nu+1:end);
p.sw.bcper=1; p=oosetfemops(p); % reassemble matrices (for fill-trafo) 
p.u=[p.mat.drop*p.u(1:p.np*p.nc.neq); aux]; 
p.nc.ilam=[3;2];      % add speed as second parameter
p.nc.nq=1;            % set auxiliary equation number to one
p.fuha.qf=@schnak_qf; % function handle for aux. eqn.
p.sw.qjac=1; p.fuha.qfder=@schnak_qfder; % analytical jac for aux. eqn.
p.nc.nsteps=5;
p.sol.ds=-0.05;
%% continue in cylinder geom. with velocity to find TW bif.
p.file.smod=10; p=cont(p,40); 
%% switch to trav.wave branch 
p=swibra('r','bpt2','t+',0.5); p.plot.bpcmp=p.nc.ilam(2); 
clf(2); p=cont(p,15); 
p=swibra('r','bpt2','t-',-0.05); p.plot.bpcmp=p.nc.ilam(2); 
p=cont(p,15);
%% plot BD of stripes 
figure(3);clf;cmp=4;ms=5; plotbraf('p',3,cmp,'cl','k'); 
plot([2.5,3.235],[2.5,3.235],'color','k','Linewidth',2);
plotbraf('q',3,cmp,'cl','b','ms',ms);
axis([1.4 3.25 2.5 4.55]);xlabel('\lambda'); ylabel('max(u)');
%% BD for cont.in domains length 
figure(3); clf; cmp=4; plotbraf('q2',3,cmp,'cl','k','ms',ms,'bplab',2);
xlabel('\rho'); ylabel('max(u)'); 
%% BD for TW cont, and soln plot 
figure(3); clf; cmp=2; ms=5; lms=5; 
plotbraf('t+','pt15',3,cmp,'cl','k','ms',ms); 
plotbraf('t-','pt15',3,cmp,'cl','k','ms',ms,'lab',15,'cl','r'); 
xlabel('\rho'); ylabel('speed');
plotsolf('t+','pt15',1,1,1); 
%% mesh-ada test, first with NBC, then with pBC 
p=loadp('q3','pt10'); p=meshada(p,'ngen',2,'sig',0.3);
%% test mesh-ada with pBC
p=loadp('t+','pt15'); p.fuha.e2rs=@e2rs_pbc; 
p.sw.verb=2; p=meshada(p,'ngen',1,'sig',0.3);
%% test ad-hoc element selection: refine nrt elements closest to xrp 
p=loadp('t+','pt15'); p.fuha.e2rs=@mye2rs; 
p.sw.verb=2; p.xrp=0; p.ntr=5; p=meshada(p,'ngen',1);
