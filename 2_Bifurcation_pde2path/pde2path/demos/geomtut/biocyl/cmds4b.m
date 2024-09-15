%% cont in lam1, at c0=0; 3 non-axi branches bifurcating at lam1<0, asymptoting lam1=0
close all; keep pphome; 
global p2pglob; p2pglob.edc='k'; p2pglob.cut=0; p2pglob.cb=1;
%% cont in lam1=tension, c0=0
mclf(2); p=swiparf('l00','pt0','0l0',2); p.nc.dsmax=0.2;  
p.sw.cdbb=1; p.sol.ds=-0.1; p.nc.foldtol=0.02; p.nc.tol=1e-6; p.np; p=cont(p,20); 
%% alternate ref. and cont
p=loadp('0l0','pt32','0l0r'); plotsol(p); sig=0.1; p.sol.ds=-0.05; p.nc.dsmax=0.05; p.sw.rlong=1; 
for i=1:3; p=refineX(p,sig); p=cont(p,10); end 
%% alternate ref. and cont, with retrigX
p=loadp('0l0r','pt45'); plotsol(p); sig=0.1; p.sol.ds=-0.05; p.nc.dsmax=0.05; p.sw.rlong=1; 
p.sw.bifcheck=0; 
for i=1:2; p=refineX(p,sig); p=retrigX(p); p=cont(p,10); end 
%% refufumaxA, doesn' really work 
p=loadp('0l0','pt32','0l0rA'); plotsol(p); p.nc.sig=0.3; p.sol.ds=-0.05; p.nc.dsmax=0.1; p.sw.rlong=1; 
p.maxA=0.02; p.fuha.ufu=@refufumaxA; p=cont(p,10); 
%% 1st non-axisym branch
p=swibra('0l0','bpt1','0l1',-0.1); p.nc.tol=1e-4; p=cont(p,5); 
%% switch on rotational PC, 
p=loadp('0l1','pt3','0l1q'); p.nc.nq=1; p.nc.ilam=[2 4]; p.sol.ds=0.1; p.nc.tol=1e-6; p=cont(p,10); 
%% alternate ref. and cont
p=loadp('0l1q','pt8','0l1qr'); plotsol(p); sig=0.075; p.sol.ds=0.05; p.nc.dsmax=0.05; 
for i=1:5; p=refineX(p,sig); p=cont(p,5); end 
%% 2nd non-axisym branch 
p=swibra('0l0','bpt2','0l2',0.1); pause; p.nc.dsmax=0.2; p.nc.tol=1e-4; p.sw.bifcheck=0; p=cont(p,3); 
%% switch on rotational PC, 
p=loadp('0l2','pt2','0l2q'); p.nc.nq=1; p.nc.ilam=[2 4]; p.sol.ds=0.1; p.nc.tol=1e-6; p=cont(p,5); 
%% alternate ref. and cont
p=loadp('0l2q','pt6','0l2qr'); plotsol(p); sig=0.075; p.sol.ds=0.05; p.nc.dsmax=0.05; 
for i=1:6; p=refineX(p,sig); p=cont(p,5); end 
%% 3rd non-axisym branch 
p=swibra('0l0','bpt3','0l3',0.1); pause; p.nc.dsmax=0.2; p.nc.tol=1e-4; p.sw.bifcheck=0; p=cont(p,3); 
%% switch on rotational PC, 
p=loadp('0l3','pt3','0l3q'); p.nc.nq=1; p.nc.ilam=[2 4]; p.sol.ds=0.1; p.nc.tol=1e-6; p=cont(p,3); 
%% alternate ref. and cont
p=loadp('0l3q','pt6','0l3qr'); plotsol(p); sig=0.075; p.sol.ds=0.05; p.nc.dsmax=0.05; 
for i=1:6; p=refineX(p,sig); p=cont(p,5); end 
%% branch plot, 
f=4; mclf(f); xlab='\lambda_1'; c=5; ylab='A'; 
%c=7; ylab='E'; %ax=[-1.9 1 -25 25]; 
%c=8; ylab='\delta_{mesh}'; ax=[-1.9 1 -25 25]; 
plotbra('0l0r','pt63',f,c,'cl','b','lab',58,'fms',0); 
plotbra('0l1qr',f,c,'cl','r','lab',10,'lab',[25,40]); 
plotbra('0l2qr',f,c,'cl','m','lab',[]); plotbra('0l3qr',f,c,'cl',p2pc('o1'),'lab',[]); 
grid on; xlabel(xlab); ylabel(ylab); 
%% soln plots
p2pglob.vi=[60,20]; p2pglob.edc='k'; p2pglob.cm='parula'; p2pglob.tsw=0; p2pglob.showbd=2; 
pplot('0l0r','pt70');  pause; pplot('0l1qr','pt40');  pause; 
pplot('0l2qr','pt60'); pause; pplot('0l3qr','pt30'); pause; 
p2pglob.edc='none'; pplot('0l1qr','pt25'); 
%% large lam1
p=loadp('0l0','pt0','0l0b'); p.sol.ds=0.2; p.nc.dsmax=2; p.nc.dlammax=2; 
p.sw.bifcheck=0; p=cont(p,20); 
p.sw.rlong=1; p.fuha.e2rs=@e2rsshape1; p=refineX(p,0.1); p=retrigX(p);mq=meshqdat(p); mq', p=cont(p,10); 
p=refineX(p,0.2); p=retrigX(p);mq=meshqdat(p); mq', p=cont(p,20); 
%% branch plot, large lam1, A or E 
f=3; mclf(f); xlab='\lambda_1'; c=5; ylab='A'; % c=7; ylab='E'; 
plotbra('0l0b','pt40',f,c,'cl','b','lab',[7 20 40]); xlabel(xlab); ylabel(ylab); 
%% branch plot, mesh-qual 
f=4; mclf(f); xlab='\lambda_1'; c=8; ylab='\delta_m'; plotbra('0l0b','pt40',f,c,'cl','b','lab',[7 20 40],'fp',1);  xlabel(xlab); ylabel(ylab); 
%%
p2pglob.vi=[10,20]; p2pglob.edc='none'; p2pglob.cm='parula'; p2pglob.tsw=0; p2pglob.showbd=2; 
pplot('0l0b','pt7'); pause; pplot('0l0b','pt20'); pause; pplot('0l0b','pt40');

