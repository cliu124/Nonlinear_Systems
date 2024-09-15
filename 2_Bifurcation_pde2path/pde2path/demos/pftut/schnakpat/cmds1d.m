close all; keep pphome; 
%% C1: init on small (arbitrary) 1D domain, and use spufu to plot disp rel. 
p=[];lx=1;nx=20; 
par=[3.5,-0.6,60]; p=schnakinit(p,lx,nx,par); p.nc.dsmax=0.5; 
p.fuha.ufu=@spufu; % set user function to "spectral plot ufu" 
p.sol.ds=-0.1; p=setfn(p,'dummy'); p=cont(p,20); % continue for plotting disp
%% C2: init on larger domain, with rather large sigma to have subcrit. stripes
p=[]; kc=sqrt(sqrt(2)-1); lx=5*2*pi/kc; nx=500; 
p=schnakinit(p,lx,nx,par); p=setfn(p,'h1D'); 
p=findbif(p,6); % many bif-points close to each other, use findbif 
p=cont(p,20);  % a few more steps (for later plotting)  
%% C3: stripes 1,2,3,6, and 1st snake on stripes 1
p=swibra('h1D','bpt1','1D1',0.1); p=cont(p);
p=swibra('h1D','bpt2','1D2',0.1); p=cont(p); 
p=swibra('h1D','bpt3','1D3',0.1); p=cont(p); 
p=swibra('h1D','bpt6','1D6',0.1); p=cont(p); 
p=swibra('1D1','bpt1','sn1D',0.1); p=cont(p,110); 
%% C4: plot BD 
fnr=3; figure(fnr); clf; c=0; plotbra('h1D','pt40',fnr,c,'cl','k'); 
plotbra('1D1',fnr,c,'cl','b', 'lab',30); 
plotbra('1D2',fnr,c,'cl',p2pc('b2'), 'lab',30); 
plotbra('1D3',fnr,c,'cl',p2pc('br1'),'lab',30); 
plotbra('1D6',fnr,c,'cl','m','lab',30); pause
plotbra('sn1D',fnr,c,'cl','r','lab',[20 60]); 
axis([2.6 3.5 30 43]); 
%% C5: sol plots  
plotsol('1D1','pt30',1,[1 2],'cl',{'k','b'}); pause; plotsol('1D2','pt30',1,1,1); pause; 
plotsol('1D3','pt30'); pause; plotsol('1D6','pt30'); pause 
plotsol('sn1D','pt20'); pause; plotsol('sn1D','pt60'); 
%% C6: fold continuation in sigma, Forced output at sig=-0.4 
figure(2); clf; p=spcontini('1D1','fpt1',2,'1D1f'); p.sol.ds=0.05; p.usrlam=-0.4; 
p.nc.lammin=-1; p.plot.bpcmp=1; p.sw.bifcheck=1; p=cont(p,40); 
p=spcontini('1D3','fpt1',2,'1D3f'); p.sol.ds=0.05; p.usrlam=-0.4; 
p.nc.lammin=-1; p.plot.bpcmp=1; p.sw.bifcheck=1; p=cont(p,30); 
p=spcontini('1D6','fpt1',2,'1D6f'); p.sol.ds=0.05; p.usrlam=-0.4; 
p.nc.lammin=-1; p.plot.bpcmp=1; p.sw.bifcheck=1; p=cont(p,40); 
%% C7: fold-cont BD plot 
fnr=4; figure(fnr); clf; c=1; plotbra('1D1f','pt40',fnr,c,'cl','b'); 
plotbra('1D3f','pt30',fnr,c,'cl',p2pc('br1')); 
plotbra('1D6f','pt40',fnr,c,'cl','m'); 
%% C8: exit fold-cont at sig=-0.4
p=spcontexit('1D1f','pt15','1D1a'); p.sol.ds=-1e-2; 
p.sw.spcalc=1; p.plot.bpcmp=0; p=cont(p,20); 
p=spcontexit('1D1f','pt15','1D1b'); p.sol.ds=1e-2; 
p.sw.spcalc=1; p.plot.bpcmp=0; p=cont(p,20); 
p=spcontexit('1D3f','pt12','1D3a'); p.sol.ds=-1e-2; 
p.sw.spcalc=1; p.plot.bpcmp=0; p=cont(p,20); 
p=spcontexit('1D3f','pt12','1D3b'); p.sol.ds=1e-2; 
p.sw.spcalc=1; p.plot.bpcmp=0; p=cont(p,20); 
p=spcontexit('1D6f','pt15','1D6a'); p.sol.ds=-1e-2; 
p.sw.spcalc=1; p.plot.bpcmp=0; p=cont(p,20); 
p=spcontexit('1D6f','pt15','1D6b'); p.sol.ds=1e-2; 
p.sw.spcalc=1; p.plot.bpcmp=0; p=cont(p,20); 
%% C9: plot BD of regular cont at sigma=-0.4
fnr=5; figure(fnr); clf; c=0; plotbra('h1D','pt40',fnr,c,'cl','k'); 
plotbra('1D1a','pt20',fnr,c,'cl','b', 'lab',30); 
plotbra('1D1b','pt20',fnr,c,'cl','b', 'lab',30); 
plotbra('1D3a','pt20',fnr,c,'cl',p2pc('br1'),'lab',30); 
plotbra('1D3b','pt20',fnr,c,'cl',p2pc('br1'),'lab',30); 
plotbra('1D6a','pt20',fnr,c,'cl','m','lab',30); 
plotbra('1D6b','pt20',fnr,c,'cl','m','lab',30);
axis([3 3.25 30 36.5]); 