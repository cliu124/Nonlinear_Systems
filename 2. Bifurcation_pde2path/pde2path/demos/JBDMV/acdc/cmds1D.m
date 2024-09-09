%% 1D Allen-Cahn dead core 
close all; keep pphome; 
%% init 
p=[]; par=[0.1   0.85   0.01 1      0.2 8 0.6];  
%          lam, ga=1/m, c,   cubic, pa, pk pm  (ampl., wave-nr, mean)    
dim=1; sw=1; p=acinit(p,1,200,par,dim,sw); p=setfn(p,'0'); plotsol(p); 
p.usrlam=[18 25]; p.Jdel=1e-8; % u_+=max(u,Jdel) in sGjac 
p=cont(p,30); 
%% center DC
p=swibra('0','bpt1','b1',-0.1); p.nc.dsmax=5; p=cont(p,100); 
% p=loadp('b1','pt20','b1r'); p=oomeshadac(p,'ngen',5,'sig',0.5); % meshada
%% further DC branches 
p=swibra('0','bpt2','b2',0.1); p.nc.almin=0.1; p.sw.jac=1; p=cont(p,100); % asymetric
p=swibra('0','bpt3','b3a',0.1); p=cont(p,80);  % center DC, 2 humps, unstable
p=swibra('0','bpt3','b3b',-0.1); p=cont(p,60);  % 2 DCs, center hump, partly stable 
%% BD, min 
f=3; c=9; figure(f); clf;  plotbra('0',f,c,'cl','k','lsw',0);
plotbra('b1',f,c,'cl','b'); 
axis([0 8 0 1.1]); ylabel('min(u)'); 
%% BD, L^2
f=3; c=0; figure(f); clf; plotbra('0',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','b','lab', [60]); plotbra('b2',f,c,'cl','m','lab',[25 60]); 
plotbra('b3a',f,c,'cl',p2pc('r1'),'lab',[10 40]); 
plotbra('b3b',f,c,'cl',p2pc('r2'),'lab',55); ylabel('||u||_2'); 
axis([0 40 0.1 1.8]); 
%% sol plot; 
mypsol('b1','pt20'); pause; mypsol('b1','pt60'); pause; mypsol('b2','pt25'); pause; 
mypsol('b2','pt60'); pause; mypsol('b3a','pt10'); pause; mypsol('b3a','pt40'); pause; 
mypsol('b3b','pt55'); 
%% cont in c: DCs lift up 
p=swiparf('b3b','pt40','b3b-c',3); p.sol.ds=0.01; p.nc.dsmax=0.1; p=cont(p,15); 
%% 2ndary bifs, 
p=swibra('b3b','bpt1','b3-1'); p=cont(p,50);