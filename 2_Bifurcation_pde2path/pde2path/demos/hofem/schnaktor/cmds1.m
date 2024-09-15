close all; keep pphome; 
%% Schnakenberg on (small) torus
R=12; rho=4; dir='6t1'; dsmin=0.01; dsmax=0.05; nx=25; Kfn='K1.mat'; try delete(Kfn); catch; end; 
p=[]; par=[3.21 -0.1 60 R rho 0]; lx=pi; ly=pi; sw.sym=2; psw=[1 2]; % full periodicity 
hofem=[]; hofem.sw=1; hofem.Kfn=Kfn; p=schnaktorinit(p,[lx,ly],nx,par,psw,sw,hofem); 
p.np, p.pm.resfac=1e-4; p.sol.ds=-dsmin; p.sw.verb=2; p.nc.bisecmax=5; 
p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.tol=1e-8; p=cont(p,4);
%% BP1, using cswibra with besw=0 to just compute many 'almost' kernel vectors ...
aux=[]; aux.m=8; aux.besw=0; % many bifs, use large m
p0=cswibra(dir,'bpt1',aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.11; p0.pm.resfac=1e-3; 
p0.u(p0.nu+2)=-0.1; p0.nc.bisecmax=5; 
%% use gentau to step through bifs, all 'simple', i.e., at most double due to transl. invar. 
p=gentau(p0,1,'c1'); p.sol.ds=-0.05; p.nc.dsmax=0.05; p=cont(p,4); pause 
p.nc.nq=1; p.nc.ilam=[1 6]; p=cont(p,20); % switch on PC and cont further
%%
p=gentau(p0,[0 0 1],'c2'); p.sol.ds=-0.05; p.sw.bifcheck=2; p=cont(p,5); 
p.nc.nq=1; p.nc.ilam=[1 6];  p=cont(p,15); 
%%
p=gentau(p0,[0 0 0 0 1],'c3'); p.sol.ds=-0.05; p.sw.bifcheck=2; p=cont(p,5); 
p.nc.nq=1; p.nc.ilam=[1 6];  p=cont(p,15); 
%%
p=gentau(p0,[0 0 0 0 0 0 1],'c4'); p.sol.ds=0.01; p=cont(p,2); 
p.nc.nq=1; p.nc.ilam=[1 6];p=cont(p,15); 
%% a secondary bif, connection between inner spots and pure rings; and 'bands around torus'  
p=swibra('c2','bpt2','c2-2',0.05); p.file.smod=5; p.sw.bifcheck=2; p=cont(p,20);  
%% plot BD 
dir='6t1'; fnr=3; figure(fnr); clf; pcmp=7; plotbra(dir,fnr,pcmp,'cl','k','lsw',0);
plotbra('c1','pt20',fnr,pcmp,'cl','m','lab',10,'fp',0); 
plotbra('c2','pt20',fnr,pcmp,'cl',p2pc('b2'),'lab',10,'fp',0); 
plotbra('c3','pt20',fnr,pcmp,'cl','r','lab',[10]); 
plotbra('c2-2','pt20',fnr,pcmp,'cl',p2pc('r2'),'lab',[5],'lp',10); 
plotbra('c4','pt20',fnr,pcmp,'cl',p2pc('b3'),'lab',20); 
axis([2.8 3.23 3.2 5]); ylabel('max u_1'); 
%% solution plots 
storplot('c1','pt10'); pause; storplot('c2','pt10'); pause; storplot('c3','pt10'); pause; 
storplot('c4','pt20'); pause; storplot('c3-2','pt5'); 
%% testing genuine Dx-matrix, has errors at bdry 
p=loadp('c1','pt10','du'); plotsol(p,1,1,2); 
ux=p.mat.Dx*p.u(1:p.nu); plotsolu(p,ux,6,1,1); 