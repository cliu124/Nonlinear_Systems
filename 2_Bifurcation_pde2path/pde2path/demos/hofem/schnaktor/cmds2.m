close all; keep pphome; 
%% Schnakenberg on larger (half-)torus
R=20; rho=10; dir='6t2'; dsmin=0.01; dsmax=0.05; nx=19; 
p=[]; par=[3.21 -0.1 60 R rho 0]; lx=pi/2; ly=pi; psw=2; % half-torus, periodic only in y; 
sw6=1; sw.sym=2; p=schnaktorinit(p,[lx,ly],nx,par,psw,sw,sw6); p.np 
p.pm.resfac=1e-4; p.sol.ds=-dsmin; p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir);
%%
p.nc.mu2=0.01; p=cont(p,10);
%% BP1
aux=[]; aux.m=8; aux.besw=0;  % many bifs, use large m, Bif.Eqns not used, just kernel 
p0=cswibra(dir,'bpt1',aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.pm.resfac=1e-3; p.sol.ds=0.1; 
p0.nc.tol=1e-6; 
%% use gentau to generate some branches near BP1, most are unstable
% use smod=2 (or even 1) to write more data to file for movies! 
p=gentau(p0,1,'6d1a'); p.sol.ds=0.05; p=cont(p,20); 
p=gentau(p0,[0  1],'6d2a'); p.sol.ds=0.05; p.file.smod=2; p=cont(p,20); 
p=gentau(p0,[0 0 1],'6d3a'); p.sol.ds=0.05;  p.file.smod=2; p=cont(p,20); 
p=gentau(p0,[0 0 0 1],'6d4a'); p.sol.ds=0.05; p.file.smod=2; p=cont(p,20);
%%
p=gentau(p0,-1,'6d1b'); p.sol.ds=0.05; p=cont(p,20); 
p=gentau(p0,[0  -1],'6d2b'); p.sol.ds=0.05; p.file.smod=2; p=cont(p,20); 
p=gentau(p0,[0 0 -1],'6d3b'); p.sol.ds=0.05;  p.file.smod=2; p=cont(p,20); 
p=gentau(p0,[0 0 0 -1],'6d4b'); p.sol.ds=0.05; p.file.smod=2; p=cont(p,20);
%% 4th BP gives stripes 
p0=cswibra(dir,'bpt3',aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.pm.resfac=1e-3;
%%
p=gentau(p0,[0 0 0 1],'63-4'); p.sol.ds=0.1; p=cont(p,30); 
%% 4th BP gives stripes 
p0=cswibra(dir,'bpt4',aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.pm.resfac=1e-3;
%%
p=gentau(p0,[0 0 0 0 1],'d10'); p.sol.ds=0.1; p=cont(p,30); 
%% plot BD 
dir='t2'; fnr=3; figure(fnr); clf; pcmp=7; plotbra(dir,'pt10',fnr,pcmp,'cl','k','lsw',0);
plotbra('6d7','pt50',fnr,pcmp,'cl','r','lab',41); %[10,41,50]); 
plotbra('6d1','pt20',fnr,pcmp,'cl','b','lab',15); 
%plotbra('d2','pt30',fnr,pcmp,'cl',p2pc('o1'),'lab',20);
%plotbra('d4','pt20',fnr,pcmp,'cl',p2pc('o3'),'lab',20);
plotbra('6d10','pt15',fnr,pcmp,'cl','k','lsw',0); %'lab',9);
axis([3 3.24 3.15 4.6]); ylabel('max u_1'); 
%% plot solns 
storplot('6d1a','pt10'); pause; storplot('6d2s','pt10',1); pause; storplot('6d4','pt10',1); pause
storplot('6d1','pt15'); pause; storplot('d10','pt9',1); pause; 
%%
storplot('d2','pt20'); pause; storplot('d4','pt20',1);
%%
p=loadp('d7','pt20'); p.file.smod=2; p=pmcont(p,50); 
%%
storplot('d1','pt15'); 
