close all; keep pphome; 
%% Schnakenberg on (small) torus
R=12; rho=4; dir='t1'; dsmin=0.01; dsmax=0.05; nx=50; 
p=[]; par=[3.21 -0.1 60 R rho 0]; lx=pi; ly=pi; sw.sym=2; psw=[1 2]; % full periodicity 
p=schnaktorinit(p,[lx,ly],nx,par,psw,sw); p.np, p.pm.resfac=1e-4; p.sol.ds=-dsmin; p.sw.verb=2; 
p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.tol=1e-6; p=cont(p,10);
%% BP1, using cswibra with besw=0 to just compute many 'almost' kernel vectors ...
aux=[]; aux.m=8; aux.besw=0; % many bifs, use large m
p0=cswibra(dir,'bpt1',aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.11; p0.pm.resfac=1e-3; 
%% use gentau to step through bifs, all 'simple', i.e., at most double due to transl. invar. 
p=gentau(p0,1,'c1'); p.sol.ds=0.05; p.nc.dsmax=0.05; p=cont(p,1); 
p.nc.nq=1; p.nc.ilam=[1 6]; p=cont(p,20); % switch on PC and cont further
p=gentau(p0,[0 0 1],'c2'); p.sol.ds=0.05; p=cont(p,1); p.nc.nq=1; p.nc.ilam=[1 6]; p=cont(p,20); 
%% c3, becomes stable, check 2ndary bifs
p=gentau(p0,[0 0 0 0 1],'c3'); p.sol.ds=0.01; p.sw.bifcheck=2; p=cont(p,1); 
p.nc.nq=1; p.nc.ilam=[1 6]; p=cont(p,20);
%% c4, unstable again 
p=gentau(p0,[0 0 0 0 0 0 1],'c4'); p.sol.ds=0.01; p=cont(p,1); p.nc.nq=1; p.nc.ilam=[1 6]; p=cont(p,20); 
%% a secondary bif, connection between inner spots and pure rings; and 'bands around torus'  
p=swibra('c3','bpt5','c3-5',0.05); p.file.smod=5; p.sw.bifcheck=0; p=cont(p,20);  
%% plot BD 
dir='t1'; fnr=3; figure(fnr); clf; pcmp=7; plotbraf(dir,'pt30',fnr,pcmp,'cl','k','lsw',0);
plotbra('c1','pt20',fnr,pcmp,'cl',p2pc('b1'),'lab',13,'fp',5); 
plotbra('c2','pt20',fnr,pcmp,'cl',p2pc('b2'),'lab',10,'fp',3); 
plotbra('c3','pt26',fnr,pcmp,'cl','r','lab',[20]); 
plotbra('c3-5','pt20',fnr,pcmp,'cl',p2pc('r2'),'lab',[8,15],'lp',20); 
plotbra('c4','pt20',fnr,pcmp,'cl',p2pc('b3'),'lab',13,'fp',2); 
axis([3.01 3.23 3.01 4.8]); ylabel('max u_1'); 
%% solution plots 
storplot('c1','pt13'); pause; storplot('c2','pt10'); pause; 
storplot('c3','pt20'); pause; storplot('c3-5','pt8'); pause; 
storplot('c3-5','pt15'); pause; storplot('c4','pt13'); 
%% continue in R: stable spots become unstable 
p=swiparf('c3','pt20','c3R',[4,1]); p.sol.ds=0.1; p.nc.dsmax=1.1; p=cont(p,20); 
%%
storplot('c3R','pt10'); 