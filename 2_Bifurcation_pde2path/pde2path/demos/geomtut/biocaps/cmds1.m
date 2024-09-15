close all; keep pphome; clear classes
%% Helfrich cap, c0=0.5; al=1; b=-0.5 (initial), then cont in c0 at different b
% see cmds1b.m for b=-4, and cmds1plot.m for plotting
global p2pglob; p2pglob.edc='k'; p2pglob.cut=0; p2pglob.showbd=2; 
p2pglob.showN=2; p2pglob.cb=0; p2pglob.vi=[20,10]; p2pglob.faceal=1; 
%% init at hemisphere 
sym=0; p=[]; srot=0;  al=1; c0=0.5; l1=0.25; b=-0.499; 
par=[al; l1; c0; b; srot]; sw=5; X=[]; tri=[]; 
p=capinit(p,sw,par); p=setfn(p,'b'); pplot(p,1); 
p.nc.dsmax=0.2; p.nc.eigref=-1000; p.nc.neig=6; p.nc.mu1=150; p.nc.mu2=0.2; 
p.nc.foldtol=0.01; p.sw.foldcheck=0; p.nc.bisecmax=8;  p.sw.cdbb=1; 
%% refine at bdry 
p.fuha.e2rs=@e2rsbdry; p.sw.nobdref=0; p.sw.rlong=1; p.tau=p.u(1:p.nu+1); 
sigr=0.2; p=refineX(p,sigr); p=retrigX(p); p0=p; 
%% axisym branch 
p=p0; p.sw.Xcont=1; p.sw.para=1; p.sol.ds=-0.01; p.nc.dsmax=0.5; 
p.nc.tol=5e-3; p=cont(p,2); 
p.nc.tol=1e-8; p=resetc(p);  p=cont(p,21); 
%% other direction
p=loadp('b','pt3','bb'); p.sol.ds=-p.sol.ds; p.nc.dsmax=0.2; p=cont(p,3); 
p=resetc(p);  p=cont(p,20); 
%% **************  cont in c0 at b=-1.66, two directions *************
mclf(2); p=swiparf('bb','pt9','c00b',3); getaux(p)', 
p.sol.ds=-0.1; p.nc.dsmax=0.2; p.sw.bifcheck=0; p=cont(p,30);
p=swiparf('bb','pt9','c00',3); p.sol.ds=0.1; p.nc.dsmax=0.2; p=cont(p,10);
%% ************** cont in c0 at b=-3.42      *************************
mclf(2); p=swiparf('bb','pt20','c01b',3); getaux(p)', p.sw.Xcont=2; pause 
p.sol.ds=-0.1; p.nc.dsmax=0.2; p=cont(p,55); pause 
p=swiparf('bb','pt20','c01',3); p.sol.ds=0.1; p.nc.dsmax=0.2; p=cont(p,10);
%% 1st non-axi 
p=swibra('c01b','bpt1','c01b-1',0.1); pplot(p); pause, 
p.nc.dsmax=0.2; p.nc.tol=1e-3; p.sw.bifcheck=0; p=cont(p,4); 
%% switch on PC, and mesh-ref
p=loadp('c01b-1','pt2','c01b-1q'); p.tau=[p.tau; 0]; 
p.nc.nq=1; p.nc.ilam=[3 5]; p.fuha.qf=@qfrot; p.fuha.qfder=@qrotder; 
p.sw.qjac=1;p.nc.tol=1e-6; p.nc.dsmax=0.2; p.sw.rlong=1; p.sw.nodref=0; 
p.fuha.e2rs=@e2rsA; p.sw.ips=2; p=cont(p,8); 
for i=1:2; p=refineX(p,0.2); p=retrigX(p); pplot(p,10); pause; p=cont(p,2*i); end 
p=cont(p,2); 
%% further refinement and cont (reload for trial and error!) 
p=loadp('c01b-1q','pt16'); 
for i=1:2; p=refineX(p,0.1); p=retrigX(p); pplot(p,10); p=cont(p,4); end 
%% 2nd non-axi, connects to 3rd BP 
p=swibra('c01b','bpt2','c01b-2',0.1); pplot(p); 
p.nc.dsmax=0.2; p.nc.tol=1e-4; p.sw.bifcheck=0; p=cont(p,4); 
%% switch on PC, 
p=loadp('c01b-2','pt4','c01b-2q'); p.tau=[p.tau; 0]; 
p.nc.nq=1; p.nc.ilam=[3 5]; p.fuha.qf=@qfrot; p.fuha.qfder=@qrotder; p.nc.eigref=-1500;
p.sw.qjac=1;p.nc.tol=1e-6; p.nc.dsmax=0.2; p.sw.bifcheck=0; p=cont(p,7); 
%% 3rd non-axi, connects to c01b-2 
p=swibra('c01b','bpt3','c01b-3',0.1); pplot(p); 
p.nc.dsmax=0.2; p.nc.tol=1e-4; p.sw.bifcheck=0; p=cont(p,9); 
%% ****************   cont in c0 at b=-4 *******************************
mclf(2); p=swiparf('bb','pt24','c02b',3); getaux(p)', p.sw.Xcont=2; pause 
p.sol.ds=-0.1; p.nc.dsmax=0.2; p=cont(p,55);
p=swiparf('bb','pt24','c02',3); p.sol.ds=0.1; p.nc.dsmax=0.2; p=cont(p,10);
%% 1st non-axi 
p=swibra('c02b','bpt1','c02b-1',0.1); pplot(p); pause, 
p.nc.dsmax=0.2; p.nc.tol=1e-3; p.sw.bifcheck=0; p=cont(p,4); 
%% switch on PC, and mesh-ref
p=loadp('c02b-1','pt4','c02b-1q'); p.tau=[p.tau; 0]; 
p.nc.nq=1; p.nc.ilam=[3 5]; p.fuha.qf=@qfrot; p.fuha.qfder=@qrotder; 
p.sw.qjac=1;p.nc.tol=1e-6; p.nc.dsmax=0.2; p.sw.rlong=1; p.sw.nodref=0; 
p.fuha.e2rs=@e2rsA; p.sw.ips=2; p=cont(p,8); 
for i=1:2; p=refineX(p,0.2); p=retrigX(p); p=cont(p,2*i); end 
%% further refinement and cont (reload for trial and error!) 
p=loadp('c02b-1q','pt18'); 
p=refineX(p,0.1); p=retrigX(p); p=cont(p,7); 
p=refineX(p,0.1); p=retrigX(p); p=cont(p,4); 
%% 2nd non-axi 
p=swibra('c02b','bpt2','c02b-2',0.1); pplot(p); 
p.nc.dsmax=0.2; p.nc.tol=1e-4; p.sw.bifcheck=0; p=cont(p,4); 
%% switch on PC, 
p=loadp('c02b-2','pt4','c02b-2q'); p.tau=[p.tau; 0]; 
p.nc.nq=1; p.nc.ilam=[3 5]; p.fuha.qf=@qfrot; p.fuha.qfder=@qrotder; p.nc.eigref=-1500;
p.sw.qjac=1;p.nc.tol=1e-6; p.nc.dsmax=0.2; p.sw.bifcheck=2; p=cont(p,8); 
%% mesh ref; reload for trial and error! 
p=loadp('c02b-2q','pt16','rt'); p.fuha.e2rs=@e2rsA; p.sw.ips=2;  p.sw.bifcheck=0; 
for i=1:2; p=refineX(p,0.2); p=retrigX(p); pplot(p,10); pause; p=cont(p,6); pause; end 
p=cont(p,2); 
%% secondary
p=swibra('c02b-2q','bpt1','c02b-2q-1',0.1); pplot(p); 
p.nc.dsmax=0.2; p.nc.tol=1e-3; p.sw.bifcheck=0; p=cont(p,10); 
%% mesh ref; reload for trial and error! 
p=loadp('c02b-2q-1','pt8','rt'); p.fuha.e2rs=@e2rsA; p.sw.ips=2;  p.sw.bifcheck=0; p.nc.tol=1e-6; 
p=refineX(p,0.2); p=retrigX(p); p=cont(p,12); 