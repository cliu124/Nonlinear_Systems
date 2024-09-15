close all; keep pphome; 
%% C1: brusselator, 1D, with per.forcing (amplitude sqrt(del) for del>0), init
p=[]; lx=3*pi/0.7; nx=100; a=3; b=9.5; du=6; dv=10; om=0.42; del=-0.1; 
par=[a b du dv del om]; p=bruinit(p,lx,nx,par); p.sw.verb=0; 
%% C2: increase del to find primary Hopf 
p=setfn(p,'hom1D'); p.nc.ilam=5; p.sol.ds=0.1; p=cont(p,3);  
%% C3: follow Hopf bif in oscillator (which also forces u) 
para=4; ds=0.1; figure(2); clf; aux.tl=30; aux.dlam=0; 
p=hoswibra('hom1D','hpt1',0.1,4,'1dh1',aux); p.hopf.bisec=5; p.nc.dsmax=0.15; 
p=setbel(p,2,1e-4,5,@lss); p.nc.tol=1e-6;  p.hopf.fltol=1e-3; p=cont(p,50); 
%% C4: PD bifs to subcritical osc. Turing pattern 
aux=[]; aux.sw=-1; p=poswibra('1dh1','bpt1','pd1',0.5,aux); p=cont(p,40); 
p=poswibra('1dh1','bpt2','pd2',0.5,aux); p=cont(p,30); 
p=poswibra('1dh1','bpt3','pd3',0.5,aux); p=cont(p,30); 
%% BD 
fnr=3; figure(fnr); clf; c=8; plotbra('1dh1','pt50',fnr,c,'cl','b'); 
plotbra('pd1','pt40',fnr,c,'cl','r','lab',20); 
plotbra('pd2','pt50',fnr,c,'cl',p2pc('r1'),'lab',20); 
plotbra('pd3','pt30',fnr,c,'cl',p2pc('r2'),'lab',20); 
xlabel('\delta'); ylabel('max(|u|)'); axis([0 0.25 4.35 7.5]); 
%% sol plots  
dirlist={'pd1','pd2','pd3'}; pt='pt20'; 
for i=1:3
  dir=dirlist{i}; hoplotf(dir,pt,1,1); v=[0 90]; 
  figure(1); title(['u1 at ' dir '/' pt]); view(v); colorbar; 
  figure(6); title(['u1(-lx,t) at ' dir '/' pt]);
  figure(8); title(['w(t) at ' dir '/' pt]); pause
end
%% Jaccheck
p.nc.del=1e-4; [Gu,Gn]=jaccheck(p);