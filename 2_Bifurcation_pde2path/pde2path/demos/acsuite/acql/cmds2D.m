%% 2D init 
p=[]; lx=1; nx=30; dim=2; par=[0.25; 0.5; 1; -0.3; 0.2]; % c0, lam, ga, del, epsi  
p=acinit(p,lx,nx,par,dim);  p=setfn(p,'2D'); p=findbif(p,2); 
%% first 2 BPs 
p=swibra('2D','bpt1','2D1a',0.1); p=cont(p,20);
p=swibra('2D','bpt1','2D1b',-0.1); p=cont(p,30); 
p=swibra('2D','bpt2','2D2a',0.1); p=cont(p,20); 
%% BD plot
fnr=3; figure(fnr); clf; pcmp=0; 
plotbra('2D',fnr,pcmp,'cl','b','lsw',0); 
plotbra('2D1a',fnr,pcmp,'cl','k','lab',[3 9 33]); 
plotbra('2D1b',fnr,pcmp,'cl','k','lsw',0); 
plotbra('2D2a',fnr,pcmp,'cl','r','lab',15); 
xlabel('\lambda'); 
%% soln plot 
plotsol('2D1a','pt3');pause; plotsol('2D1a','pt9');pause; 
plotsol('2D1a','pt30');pause; plotsol('2D2a','pt15');
%% a mesh-adaption test
p=loadp('2D2a','pt15'); plotsol(p); pause; p.nc.ngen=1; p=oomeshada(p); plotsol(p,6,1,1); 
%% jaccheck; rel error is about 0.001, but cont works  
[Gu,Gn]=jaccheck(p); Gd=abs(Gu-Gn); e1=max(max(Gd)); figure(10); spy(Gd>e1/2); 