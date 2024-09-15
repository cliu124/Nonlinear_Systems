close all; keep pphome; 
global p2pglob; p2pglob.edc='k'; p2pglob.cut=0;  p2pglob.tsw=0;
%% Helfrich cylinder, continuation in al 
al=1; lx=1; ly=pi; nx=35; ny=round(2*al*nx);  sym=0; p=[]; l1=10; c0=0; srot=0; 
par=[al; l1; c0; srot]; 
p=cylinit(p,lx,nx,ny,par,sym); p=setfn(p,'ab'); pplot(p,1); p.sw.Xcont=2; p.sw.cdbb=1; 
p.nc.dsmax=0.2; p.nc.eigref=-50; p.nc.mu1=5; p.nc.mu2=0.2; p.nc.bisecmax=8;  
%% decreasing al, 
p.nc.ilam=1; p.sol.ds=-0.1; p=cont(p,10); 
%% refine, retrig, cont on c0 
p=loadp('ab','pt6','abr'); plotsol(p); sigr=0.1; mq=meshqdat(p); mq' 
p.fuha.e2rs=@e2rsshape1; p.sw.nobdref=0; p.sw.rlong=1; p.DIR=[]; int=3;
for i=1:int;
p=refineX(p,sigr); mq=meshqdat(p); mq', p=retrigX(p); mq=meshqdat(p); mq' 
p.nc.dsmax=0.1; p=cont(p,5); 
end 
%% other direction, 
p=loadp('ab','pt0','a'); p.sol.ds=-p.sol.ds; p.nc.dsmax=2; p=cont(p,10); 
%% branch plot, A and E 
f=3; mclf(f); xlab='\alpha'; c=5; ylab='A'; c=7; ylab='E';  %c=9; ylab='\delta'; 
plotbra('a',f,c,'cl','b','lab',10,'fp',0,'lp',51);  
plotbra('abr',f,c,'cl',p2pc('b2'),'lab',[8 40],'fp',0,'lp',40);  
xlabel(xlab); ylabel(ylab); grid on
%%
p2pglob.cut=4; p2pglob.cm='parula'; p2pglob.vi=[10,20]; p2pglob.edc='k';
pplot('abr','pt8'); pause; pplot('abr','pt40'); pause; 
p2pglob.cut=0;p2pglob.edc='none'; pplot('a','pt10'); 


