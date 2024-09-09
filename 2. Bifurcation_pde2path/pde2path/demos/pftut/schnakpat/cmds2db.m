%% commands for Schnakenberg on longer domain, with snaking
close all; keep pphome; 
p=[]; kc=sqrt(sqrt(2)-1); lx=8*pi/kc; ly=lx/sqrt(3)/4; nx=150; par=[3.3, 0, 60]; 
sw.sym=2; p=schnakinit(p,[lx,ly],nx,par,sw); % init with criss-cross mesh 
p.fuha.outfu=@schnakbra; p=setfn(p,'bhom'); p.file.smod=5; p=cont(p,10); 
%% hex+ via qswibra
aux=[]; aux.m=2; p0=qswibra('bhom','bpt1',aux); p0.sol.ds=-0.1; p0.sw.bifcheck=1;
p0.nc.lammin=2.5; p=seltau(p0,3,'bh+',2); p=pmcont(p,30); 
%% stripes+ via gentau 
p=gentau(p0,1); p=setfn(p,'bs+'); p.nc.dsmax=0.05; p.sw.bifcheck=2; p=pmcont(p,40);
%% b+ 
p=swibra('bs+','bpt3','bb+',-0.05); p.nc.dsmin=0.05; p.nc.dsmax=0.1; p.sw.bifcheck=2;  
p=findbif(p,2); p=cont(p,15);
%% snake of loc. hex over stripes 
ds=0.02; p=swibra('bb+','bpt2','sn1',-ds); p.nc.dsmin=ds; p.nc.dsmax=2*ds;  
p.sw.bifcheck=0; p=pmcont(p,80);
%% plot BD 
fnr=3; figure(fnr); clf; cmp=6; 
plotbra('bhom','pt10',fnr,cmp,'cl','k'); plotbra('bs+',fnr,cmp,'cl','b'); 
plotbra('bh+',fnr,cmp, 'cl','m'); plotbra('bb+',fnr,cmp,'cl','r','lab',10,'lp',20); 
plotbra('sn1','pt50',fnr,cmp, 'cl','k','lab',[25 35],'lp',50); 
axis([2.5 3.25 3.18 4.1]); xlabel('\lambda'); ylabel('||u||_8'); 
%% soln plot 
p=loadp('bb+','pt5'); plotsol(p); colormap hot; fouplot(p,10,1,[20,3],'sn1/pt10'); pause
p=loadp('bb+','pt10'); plotsol(p); colormap hot; fouplot(p,10,1,[20,3],'bb+/pt10'); 