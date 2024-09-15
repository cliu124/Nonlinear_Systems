close all; keep pphome; 
global p2pglob; p2pglob.edc='k'; p2pglob.cut=0;  p2pglob.tsw=0;
%% Helfrich cylinder, continuation in c0, (alpha,lam_1)=(1,1/4)
al=1; lx=1; ly=pi; nx=30; ny=round(2*al*nx);  sym=0; p=[]; l1=10; c0=0; srot=0; 
par=[al; l1; c0; srot]; 
p=cylinit(p,lx,nx,ny,par,sym); p=setfn(p,'c0b'); pplot(p,1); p.sw.Xcont=2; p.sw.cdbb=1; 
p.nc.dsmax=0.2; p.nc.eigref=-50; p.nc.mu1=5; p.nc.mu2=0.2; p.nc.bisecmax=8;  
%% decreasing c0, with refinement, e2rs=@e2rsshape1, inradius/max(edge-lengths) 
p.fuha.e2rs=@e2rsshape1; p.sw.nobdref=1; p.sw.rlong=1; p.DIR=[]; p.nc.dsmax=0.4; 
p.sol.ds=-0.1; p.nc.tol=1e-6; tic; p=cont(p,25); toc
sigr=0.2; p=refineX(p,sigr); mq=meshqdat(p); mq', p=retrigX(p); mq=meshqdat(p); 
p=cont(p,10); mq=meshqdat(p); mq', 
%% refine again and cont further 
sigr=0.2; p=refineX(p,sigr); mq=meshqdat(p); mq', p=retrigX(p); mq=meshqdat(p); mq'
p=cont(p,5); mq=meshqdat(p); mq', 
%% other direction, 
p=loadp('c0b','pt0','c0'); p.sol.ds=-p.sol.ds; p.nc.dsmax=0.1; p=cont(p,25); 
%% refine, retrig, cont on c0 
p=loadp('c0','pt25','c0r'); plotsol(p); sigr=0.2; mq=meshqdat(p); mq' 
p.fuha.e2rs=@e2rsshape1; p.sw.nobdref=1; p.sw.rlong=1; p.DIR=[]; 
p=refineX(p,sigr); mq=meshqdat(p); mq', p=retrigX(p); mq=meshqdat(p); mq' 
p.nc.dsmax=0.1; p=cont(p,8); 
%% twice more 
sigr=0.15; p.nc.dsmax=0.05; p.sw.bifcheck=0; 
p=refineX(p,sigr); mq=meshqdat(p); mq', p=retrigX(p); mq=meshqdat(p); mq', 
p=cont(p,12); p=refineX(p,sigr); mq=meshqdat(p); mq',  p=retrigX(p); mq=meshqdat(p); mq', 
p=cont(p,5); 
%% 1st non-axisym branch 
p=swibra('c0','bpt1','c1',0.1); pause; p.nc.dsmax=0.1; p.nc.tol=1e-6; pause 
p.sw.para=2; p.sw.bifcheck=0; p=cont(p,4); 
%% refine, retrig, cont on c1
p=loadp('c1','pt1','c1r'); plotsol(p); sigr=0.1; mq=meshqdat(p); mq' 
p.fuha.e2rs=@e2rsshape1; p.sw.nobdref=1; p.sw.rlong=1; p.DIR=[]; 
p=refineX(p,sigr); mq=meshqdat(p); mq', p=retrigX(p); mq=meshqdat(p); mq', pplot(p); 
p.nc.dsmax=0.1; p=cont(p,5); 
%% branch plot, A 
f=3; mclf(f); xlab='c_0'; c=5; ylab='A'; %c=7; ylab='E'; 
plotbra('c0br','pt45',f,c,'cl',p2pc('b2'),'lab',[20 30 45],'fp',1);  
xlabel(xlab); ylabel(ylab); 
%% branch plot, mesh-qual 
f=4; mclf(f); xlab='c_0'; c=8; ylab='\delta'; 
plotbra('c0br','pt45',f,c,'cl',p2pc('b2'),'lab',[25 35 45],'fp',1);  
xlabel(xlab); ylabel(ylab); 
%% branch plot, mesh-qual, c0>0 
f=3; mclf(f); xlab='c_0'; c=8; ylab='\delta'; 
plotbra('c0r','pt47',f,c,'cl','b','lab',[33 46],'fp',0); 
plotbra('c1r','pt7',f,c,'cl','r','lab',[5]);
xlabel(xlab); ylabel(ylab); 
%% branch plot, c0>0, A or E 
f=4; mclf(f); xlab='c_0'; c=5; ylab='A';  c=7; ylab='E'; 
plotbra('c0r','pt47',f,c,'cl','b','lab',[10 25 33 46],'fp',0); 
plotbra('c1r','pt7',f,c,'cl','r','lab',[5]);
xlabel(xlab); ylabel(ylab); 
%% soln plots
p2pglob.vi=[10,20]; p2pglob.edc='k'; p2pglob.cm='spring'; p2pglob.showbd=2; p2pglob.cb=0; 
plotsol('c0','pt10'); pause; plotsol('c0','pt25'); pause;
plotsol('c0r','pt26');pause;plotsol('c0r','pt35');pause; plotsol('c0r','pt47'); pause 
plotsol('c1r','pt5');pause; 
%%
p2pglob.edc='k';  
plotsol('c1r','pt5'); pause  
plotsol('c0','pt25'); pause; plotsol('c0r','pt35');pause; plotsol('c0r','pt41'); pause 
%%
plotsol('c0b','pt20');  pause; plotsol('c0br','pt30'); pause; 
plotsol('c0br','pt35'); pause; plotsol('c0br','pt36'); pause; 
plotsol('c0br','pt45'); 