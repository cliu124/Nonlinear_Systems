%% Schnakenberg on Watts-Strogatz graph 
close all; keep pphome; 
%% init 
p=[]; par=[3.5,-0.3,1,500]; 
sw=2; np=100;   sw=0; np='G1'; % uncomment to use previously saved graph 
p=schnakinitG(p,par,np,sw); p=setfn(p,'tr1'); p.nc.lammax=12; p.nc.dsmax=0.25; 
p.sw.bifcheck=2; p.sw.verb=2; p.nc.mu1=2; p.nc.mu2=0.5; p.sol.ds=-0.1; p.nc.ilam=1; 
p.fuha.outfu=@hobra; p.plot.bpcmp=6; p.vol=1; pause; p=cont(p,25);
%% swibra   
p=swibra('tr1','bpt1','b1-1',0.1); pause; p.nc.lammin=-10; p=cont(p,50);  
%% hoswibra 
para=4; ds=0.2; tl=50; figure(2); clf;  aux=[]; aux.tl=tl; aux.dlam=0; 
p=hoswibra('b1-1','hpt1',ds,para,'h1',aux); nsteps=15; flc=1; 
p.hopf.jac=1; p.nc.dsmax=0.5; p.hopf.xi=0.05; p.file.smod=2; p.sw.verb=2; 
p.hopf.flcheck=flc; p.sw.bifcheck=1;  p.nc.tol=1e-6; 
p.hopf.pind=[1 round(tl/2)]; p.hopf.lay=[1 2]; 
t1=tic; p=cont(p,nsteps); toc(t1) 
%% reload, set smaller step size for better PD localization 
p=loadp('h1','pt14'); p.nc.dsmax=0.1; p.sol.ds=0.1; p=cont(p,10);
%% cont again after PD with larger ds 
p=loadp('h1','pt24'); p.nc.dsmax=0.5; p.sol.ds=0.5; p=cont(p,10);
%% PD 
huclean(p); ds=0.2; p=poswibra('h1','bpt1','pd1',ds); p.nc.tol=1e-4; p.nc.dsmax=2; 
p.sw.bifcheck=0; p.hopf.flcheck=1; p.nc.dsmin=1e-3; p.sw.verb=0; p.hopf.fltol=1e-2; 
p=cont(p,10); %p.nc.tol=1e-2; p=cont(p,5); p.nc.tol=1e-6; p=cont(p,23); 
%% BD plot 
f=3; c=6; figure(f); clf; plotbra('tr1',f,c,'cl','k','lsw',0); 
plotbra('b1-1',f,c,'cl','b','lab',10); plotbra('h1','pt36',f,c,'cl','r','lab',36); 
plotbra('pd1','pt10',f,c,'cl',p2pc('o1'),'lab',10); ylabel('max u'); 
%% solution plots 
plotsol('b1-1','pt10'); pause; plotsol('b1-1','pt20'); pause; 
p=loadp('h1','pt36'); incr=11; 
p.hopf.pind=1:incr:4*incr-1; p.hopf.lay=[2 2]; hoplot(p,1,1); pause 
p=loadp('pd1','pt10'); incr=18;  %incr=32; 
p.hopf.pind=1:incr:4*incr-1; p.hopf.pind=[1 35 55 80]; 
p.hopf.lay=[2 2]; hoplot(p,1,1);