% driver script for CPs for vegOC, set globals 
global s0 s1 u0 u1 Psi par xi um1 um2 sig;  
%% R=26, CP from p1/p11 to FSS, set some options, then call iscnat and plot
opt=[]; opt=ocstanopt(opt); opt.rhoi=1; opt.t1=200; opt.start=1; 
opt.tv=[]; opt.nti=20; opt.Nmax=200; opt.msw=1; opt.Stats_step='on';  
sd0='FSS'; sp0='pt18'; sd1='p1'; sp1='pt11'; flip=1; % from p1 to FSS (R=26)
fn=setfnflip(sd0,sp0,sd1,sp1,flip); 
%%
alvin=[0.1 0.25 0.5 0.75 1]; v=[15,30];
[alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,[],[],opt,fn); 
figure(15); clf; vegsolplot(sol,v,4,fn); 
%% R=26, CP from p1/p11 to p1/pt38,  
opt=[]; opt=ocstanopt(opt); opt.rhoi=1; opt.t1=200; opt.start=1; 
opt.tv=[]; opt.nti=20; opt.Nmax=200; opt.msw=1; opt.Stats_step='on';  
sd0='p1'; sp0='pt11'; sd1='p1'; sp1='pt38'; flip=0; 
fn=setfnflip(sd0,sp0,sd1,sp1,flip); alvin=[0.1 0.25 0.5 0.75 1]; v=[15,30];pause
[alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,[],[],opt,fn); 
figure(15); clf; vegsolplot(sol,v,4,fn); 
%% R=26, CP from FSS to p1/pt38,  
opt=[]; opt=ocstanopt(opt); opt.rhoi=1; opt.t1=200; opt.start=1; 
opt.tv=[]; opt.nti=20; opt.Nmax=200; opt.msw=1; opt.Stats_step='on';  
sd0='FSS'; sp0='pt18'; sd1='p1'; sp1='pt38'; flip=0; 
fn=setfnflip(sd0,sp0,sd1,sp1,flip); alvin=[0.1 0.25 0.5 0.75 1]; v=[15,30];
opt.AbsTol=1e-2; 
[alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,[],[],opt,fn); 
figure(15); clf; vegsolplot(sol,v,4,fn); 
%% R=20: 
opt.nti=10; sd0='FSS'; sp0='pt29'; sd1='p1'; sp1='pt49'; flip=0; % from FSS to p1(high) 
fn=setfnflip(sd0,sp0,sd1,sp1,flip); 
[alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,[],[],opt,fn); 
figure(15); clf; vegsolplot(sol,v,15,fn); 
%% R=10: 
sd0='FSS'; sp0='pt45'; sd1='p1'; sp1='pt65'; flip=0; % from FSS to p1(high) 
fn=setfnflip(sd0,sp0,sd1,sp1,flip); opt.msw=1; opt.Itnlmax=10; 
[alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,[],[],opt,fn); vegsolplot(sol,v,10,fn); 
figure(15); clf; vegsolplot(sol,v,15,fn);