%% refine on bN1, good
global p2pglob; 
p2pglob.pbctol=1e-3; p=loadp('bN1','pt10','bN1r'); p.sw.Xfill=0; 
plotsol(p); p.sw.orgper=0; sig=0.1; p.sw.rlong=1; p=refineX(p,sig); pause; 
p.sol.ds=0.1*sign(p.sol.ds); p.nc.tol=1e-3; p=cont(p,10); 
%% after 10 steps refine again 
p=refineX(p,sig); p=cont(p,10); 
%% refine on bN2, good
p=loadp('bN2','pt10','bN2r'); plotsol(p); p.sw.Xfill=0; 
sig=0.1; p.sw.rlong=1; p=refineX(p,sig); pause; 
p.sol.ds=0.1*sign(p.sol.ds); p.nc.tol=0.001; p=cont(p,10); 
%% moveX; better mesh
p=loadp('bN1r','pt20','bN1rm'); plotsol(p); 
delt=0.01; nit=2; p=moveX(p,delt,nit); pause; 
p.sol.ds=0.1*sign(p.sol.ds); p.nc.tol=0.001; p=cont(p,10); 
%% retrigX; same as move 
p=loadp('bN1r','pt20','bN1rrt'); pplot(p,1); p=retrigX(p); pplot(p,10); pause; 
p.sol.ds=0.1*sign(p.sol.ds); p.nc.tol=0.001; p=cont(p,10); 
%% coarsen bN2; works OK 
p=loadp('bN2','pt10','bN2c'); plotsol(p); p2pglob.pbctol=5e-3; 
keepbd=1; sig=0.5; nit=6; p=degcoarsenX(p,sig,nit,keepbd); pause; 
p.sol.ds=0.1*p.sol.ds; % p.sol.ds=-0.1*sign(p.sol.ds); 
p.nc.tol=1e-4; p=cont(p,10); 
%% refine again, works 
p=loadp('bN2c','pt20','bN2cr'); plotsol(p); p2pglob.pbctol=5e-3; p.sw.nobdref=1; 
sig=0.1; p.sw.rlong=1; p=refineX(p,sig); pause; p=cont(p,5); 
p=refineX(p,sig); pause; p=cont(p,5);  % and once more 
%% branch plot
f=3; mclf(f); p2pglob.pbctol=1e-4;  ylab='r'; xlab='\delta'; 
xlab='\delta';  c=[6 19]; %xlab='z'; ylab='A'; 
plotbra('bN1rrt','pt30',f,c,'cl','r','labi',[10]); 
plotbra('bN2cr','pt30',f,c,'cl','b','labi',[10]); 
xlabel(xlab); ylabel(ylab);  
%% 
p2pglob.cut=0; plotsol('bN1','pt10'); pause; plotsol('bN1rrt','pt25');
%%
p2pglob.cut=1; plotsol('bN2','pt20'); pause;
plotsol('bN2c','pt20'); pause;  plotsol('bN2cr','pt30'); 
