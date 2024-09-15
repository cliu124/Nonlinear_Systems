close all; keep pphome; 
global p2pglob; p2pglob.edc='k'; p2pglob.cut=0;  p2pglob.tsw=0;
%% Helfrich cylinder, continuation in c0 at larger/smaller alpha 
al=1.4; lx=1; ly=pi; nx=30; ny=round(2*al*nx);  sym=0; p=[]; l1=10; c0=0; srot=0; 
par=[al; l1; c0; srot]; p=cylinit(p,lx,nx,ny,par,sym); p=setfn(p,'b/c0b'); p.sw.cdbb=1; 
p.nc.dsmax=0.2; p.nc.eigref=-50; p.nc.mu1=5; p.nc.mu2=0.2; p.nc.bisecmax=8;  
%% decreasing c0, with refinement, e2rs=@e2rsshape1, inradius/max(edge-lengths) 
p.fuha.e2rs=@e2rsshape1; p.sw.nobdref=1; p.sw.rlong=1; p.DIR=[]; p.nc.dsmax=0.4; 
p.sol.ds=-0.1; p.nc.tol=1e-6; tic; p=cont(p,1); toc, p.up=p.u; 
sigr=0.2; p=refineX(p,sigr); p=retrigX(p); mq=meshqdat(p); p=cont(p,10); mq=meshqdat(p); mq', 
%% refine again and cont further 
sigr=0.2; p=refineX(p,sigr); p=retrigX(p); mq=meshqdat(p); mq', 
p=cont(p,5); mq=meshqdat(p); mq', 
%% other direction, 
p=loadp('b/c0b','pt0','b/c0'); p.sol.ds=-p.sol.ds; p.nc.dsmax=0.1; p=cont(p,35); 
%% refine, retrig, cont on c0 
p=loadp('b/c0','pt35','b/c0r'); plotsol(p); sigr=0.1; mq=meshqdat(p); mq' 
p.fuha.e2rs=@e2rsshape1; p.sw.nobdref=1; p.sw.rlong=1; p.DIR=[]; int=3;
for i=1:int;
  p=refineX(p,sigr); p=retrigX(p); mq=meshqdat(p); mq', p.nc.dsmax=0.1; p=cont(p,5); 
end 
%% careful near fold 
p=loadp('b/c0r','pt45','b/c0r'); plotsol(p); sigr=0.1; mq=meshqdat(p); mq' 
p.fuha.e2rs=@e2rsshape1; p.sw.nobdref=1; p.sw.rlong=1; int=7; p.sol.ds=0.01; 
p.nc.dsmax=0.05; p.sw.bifcheck=0; 
%%
for i=1:int; p=refineX(p,sigr); p=retrigX(p); p=cont(p,5); end 
%% 1st non-axisym branch 
p=swibra('b/c0r','bpt1','b/c1',-0.1); p.nc.dsmax=0.1; p.nc.tol=1e-4; pause 
p.sw.para=2; p.sw.bifcheck=0; p=cont(p,4); 
%% switch on rotational PC, 
p=loadp('b/c1','pt2','b/c1q'); p.nc.nq=1; p.nc.ilam=[3 4]; p.nc.tol=2e-4; p=cont(p,6); 
%% refine 
p=loadp('b/c1q','pt4','b/c1qr'); pplot(p); pause; sigr=0.1; mq=meshqdat(p); mq' 
p.fuha.e2rs=@e2rsshape1; p.sw.nobdref=1; p.sw.rlong=1; int=10; p.nc.dsmax=0.05; 
p.sw.bifcheck=0; p.nc.tol=1e-4; 
for i=1:int; p=refineX(p,sigr); p=retrigX(p); mq=meshqdat(p); mq', p=cont(p,4); end 
%% branch plots  
f=3; mclf(f); xlab='c_0'; c=5; ylab='A'; c=7; ylab='E';% c=9; ylab='\delta'; 
plotbra('b/c0b',f,c,'cl',p2pc('b2'),'lp',8);  
plotbra('b/c0r',f,c,'cl','b','lab',[20 51],'lp',52);  
plotbra('b/c1qr',f,c,'cl',p2pc('r2'),'lab',15); 
xlabel(xlab); ylabel(ylab); 
%% soln plots 
p2pglob.vi=[10,20]; p2pglob.edc='none'; p2pglob.cm='parula'; p2pglob.showbd=2; p2pglob.cut=0; 
pplot('b/c0','pt20'); pause; pplot('b/c0r','pt51'); pause;  pplot('b/c1qr','pt15'); pause 
%% alpha=0.6, 'small' 
al=0.6; lx=1; ly=pi; nx=30; ny=round(2*al*nx);  sym=0; p=[]; l1=10; c0=0; srot=0; 
par=[al; l1; c0; srot]; p=cylinit(p,lx,nx,ny,par,sym); p=setfn(p,'s/c0'); p.sw.cdbb=1; 
p.nc.dsmax=0.2; p.nc.eigref=-50; p.nc.mu1=5; p.nc.mu2=0.2; p.nc.bisecmax=8; 
p.fuha.e2rs=@e2rsshape1; p.sw.nobdref=1; p.sw.rlong=1; p.nc.dsmax=0.4; 
p.sol.ds=0.1; p.nc.tol=1e-6; p=cont(p,10); 
%% continue with refinement 
p=loadp('s/c0','pt10','s/c0r'); plotsol(p); sigr=0.1; int=10; p.sol.ds=0.01; p.nc.dsmax=0.1; 
for i=1:int; p=refineX(p,sigr); p=retrigX(p); mq=meshqdat(p); mq', p=cont(p,5); end 
%% 1st non-axisym branch, hard due to mesh contraction 
p=swibra('s/c0r','bpt1','s/c1',-0.1); p.nc.dsmax=0.1; p.nc.tol=1e-4; pause 
p.sw.para=2; p.sw.bifcheck=0; p=cont(p,5); 
%% branch plot, negative c0 
f=3; mclf(f); xlab='c_0'; c=5; ylab='A'; c=7; ylab='E'; %c=9; ylab='\delta'; 
plotbra('s/c0b',f,c,'cl',p2pc('b2'),'lab',[10 30]);  
xlabel(xlab); ylabel(ylab); 
%% branch plots, positive c0 
f=3; mclf(f); xlab='c_0'; c=5; ylab='A'; c=7; ylab='E'; %c=9; ylab='\delta'; 
plotbra('s/c0r',f,c,'cl','b','lab',[15 20 30],'lp',32);  
plotbra('s/c1',f,c,'cl',p2pc('r2'),'lab',5); 
xlabel(xlab); ylabel(ylab); 
%%
p2pglob.vi=[10,20]; p2pglob.edc='none'; p2pglob.cm='parula'; p2pglob.showbd=2;p2pglob.cut=0; 
pplot('s/c0b','pt10'); pause; pplot('s/c0b','pt30'); pause;  pplot('s/c1','pt5'); pause 
pplot('s/c0r','pt15'); pause; pplot('s/c0r','pt20'); pause;  pplot('s/c0r','pt30'); 