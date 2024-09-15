%% ac2D cmds1; 
close all; keep pphome; 
%% first run with standard (3-node triangles) 
p=[]; par=[1 0.2 1 0]; lx=2*pi; ly=pi; nx=20; bpoints=[5/16, 1/2, 13/16]; 
hofem.sw=0; hofem.t2sinterpol=0; p=acinit(p,lx,ly,nx,par,hofem); 
p.nc.ilam=2; p.nc.lammax=1; p.sol.ds=0.1; p.nc.dsmax=0.1; p=setfn(p,'tr'); 
p=findbif(p,3);
%% switch to first 3 bifurcating branches and continue
p=swibra('tr','bpt1','b1',0.1); b(1)=getlam(p); p=cont(p); 
p=swibra('tr','bpt2','b2',0.1); b(2)=getlam(p); p=cont(p); 
p=swibra('tr','bpt3','b3',0.1); b(3)=getlam(p); p=cont(p); 
for i=1:3; fprintf('%g ',bpoints(i)-b(i)); end; % error in BPs 
fprintf(',   ||bpoints-b||_2=%g \n', norm(bpoints-b)); 
%% reload pt0, switch on 2nd order FEM, find BPs 
p=loadp('tr','pt0','trs'); p.hofem.sw=1; p.sol.ds=0.1; p=oosetfemops(p); p=findbif(p,3);
%% switch to first 3 bifurcating branches and continue
p=swibra('trs','bpt1','b1s',0.1); b6(1)=getlam(p); p=cont(p); 
p=swibra('trs','bpt2','b2s',0.1); b6(2)=getlam(p); p=cont(p); 
p=swibra('trs','bpt3','b3s',0.1); b6(3)=getlam(p); p=cont(p); 
for i=1:3; fprintf('%g ',bpoints(i)-b6(i)); end;  % error in BPs 
fprintf(',   ||bpoints-b||_2=%g \n', norm(bpoints-b6)); 
%% plot BD 
f=3; c=0; figure(f); clf; plotbra('trs',f,c,'cl','k'); 
plotbra('b1s',f,c,'cl','b','lab',10); plotbra('b2s',f,c,'cl','r','lsw',0); 
plotbra('b3s',f,c,'cl','m','lsw',0); xlabel('\lambda'); ylabel('||u||_2'); 
%%
plotsol('b1s','pt10'); 
%% compare 3-node and 6-node K 
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
figure(1); clf; spy(K); pause; figure(6); clf; spy(p.mat.K); 
%% continue in some other param, here the BC coefficent d 
p=swiparf('b1s','pt10','b1s-dc',4); p.sol.ds=0.1; p.nc.lammin=-1; p.nc.lammax=2; 
p=cont(p,15); 
%% mesh-adaption 
p=loadp('b1s-dc','pt10'); plotsol(p,1,1,1); p=setfn(p,'10r'); pause 
p.hofem.sw=0; % use 3-node triangles (same nodal vals) 
p=resetc(p); p.plot.pstyle=1; plotsol(p); p.np, 
p.nc.ngen=1; p=oomeshada(p); plotsol(p);  stansavefu(p);
%% switch back to 6-node triangles 
p.hofem.sw=1; p.hofem.t2sinterpol=1; % interpolate 3-node soln to new 6-node mesh 
p=tri2six(p); p.pdeo.grid=setidssq(p.pdeo.grid); 
p=setfemops(p); p=setfn(p,'dcr');  stansavefu(p); plotsol(p); 
%%
huclean(p); p=cont(p,10); 
