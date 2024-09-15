%% CSS  SLOC2D 
close all; p=[];lx=2*pi/0.44; ly=lx/2; nx=20; ndim=2; 
sw=1; p=slinit(p,lx,ly,nx,sw,ndim); p=setfn(p,'2DFSC'); 
%% find bif.points from FSC/FSI branch 
p.nc.nsteps=100; p=cont(p); 
%% bif from F2 (set bpt* and h* and repeat as necessary) 
p=swibra('2DFSC','bpt2','h2',0.05); p.nc.dsmax=0.3; 
p.sw.spcalc=0; p.sw.bifcheck=0; p=cont(p,80); 
p=swibra('2DFSC','bpt3','h3',0.05); p.nc.dsmax=0.3; 
p.sw.spcalc=0; p.sw.bifcheck=0; p=cont(p,80); 
%% Bifurcation Diagram L2  
figure(3); clf; pcmp=3; ms=5; 
plotbra('2DFSC',3,pcmp,'ms',ms,'cl','k'); 
plotbra('h2',3,pcmp,'ms',ms,'cl','r'); 
plotbra('h3',3,pcmp,'ms',ms,'cl','c'); 
[lh,oh,ph]=legend('FSS','h2','h3');
set(ph(2),'color','r'); set(ph(3),'color','c'); 
xlabel('b'); ylabel('||P ( )||_2'); %axis([0.6 0.75 0.4 0.95]);
%%
plotsolf('h2','pt16',1,1,2); xlabel(''); ylabel(''); axis image; pause 
plotsolfu('h2','pt16',1,2,2); xlabel(''); ylabel(''); axis image; colormap hot; pause 
plotsolf('h3','pt15',1,1,2);xlabel(''); ylabel(''); axis image; pause 
plotsolfu('h3','pt15',1,2,2);xlabel(''); ylabel(''); axis image; colormap hot;pause 
%%
cssvalf('h2','pt15'); 