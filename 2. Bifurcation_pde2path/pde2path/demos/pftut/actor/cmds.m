%% AC on torus 
close all; keep pphome; 
%% init and continue trivial branch 
p=[]; R=2; rho=1; par=[R rho -0.3 1 0]; % parameters [R r lambda gamma s]
lx=pi; ly=pi; p=acinit(p,lx,ly,30,par); p.sw.bifcheck=2; p.nc.mu1=0.5; p.nc.mu2=0.1; 
p.sw.verb=1; p.nc.bisecmax=10; p.nc.neig=30; 
p.nc.ilam=3; p.sol.ds=0.1; p.nc.dsmax=0.11; p=setfn(p,'tr'); 
tic; p=cont(p,35); toc
%% spat.hom.branch, transcrit, hence both directions 
p=swibra('tr','bpt1','b1a',0.05); p=cont(p,30); 
p=swibra('tr','bpt1','b1b',-0.05); p=cont(p,30); 
%% other simple BPS, cause no phi-dep, b5 is transcritical 
p=swibra('tr','bpt4','b4',0.1); pause; p=cont(p,30); 
p=swibra('tr','bpt5','b5a',0.1); pause; p=cont(p,30); 
p=swibra('tr','bpt5','b5b',-0.1); pause; p=cont(p,30); 
%% others (double, but some effectively simple due to symmetry, x-dependence, hence need PC)
aux=[]; aux.m=2;  
bpli=[2 3 6 7]; 
nbp=length(bpli); 
for i=1:nbp; ii=bpli(i); bp=['bpt' mat2str(ii)]; br=['b' mat2str(ii)]; 
  p=swibra('tr',bp,br,0.1); p=cont(p,4); % a few steps without PC 
  p.nc.nq=1; p.sol.ds=0.05; p.nc.ilam=[3 5]; p=cont(p,40); % set PC 
end 
%% bifurcation diagram plotting, 
f=3; c=0; figure(f); clf; plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1a',f,c,'cl',p2pc('r1'),'lsw',0); plotbra('b1b',f,c,'cl',p2pc('r1'),'lsw',0);
plotbra('b2',f,c,'cl',p2pc('b1')); plotbra('b3',f,c,'cl',p2pc('b3')); 
plotbra('b4',f,c,'cl',p2pc('r2'),'lsw',0); 
plotbra('b6',f,c,'cl',p2pc('o1')); 
plotbra('b5a',f,c,'cl',p2pc('r3')); plotbra('b5b',f,c,'cl',p2pc('r3'));
plotbra('b7',f,c,'cl',p2pc('o2')); 
axis([-0.3 2.15 0 8]); xlabel('\lambda'); ylabel('||u||_2'); 
%% solution plots 
torplot('b2','pt30'); pause; torplot('b3','pt23'); pause 
torplot('b4','pt16'); pause; torplot('b5a','pt14');  pause ; torplot('b5b','pt15');  pause 
torplot('b6','pt11'); pause; torplot('b7','pt8');  
%% continue in R 
p=swiparf('b2','pt30','b2R',[1 5]); p.sol.ds=0.1; p.nc.lammax=5.1; huclean(p); p=cont(p,40);
p=swiparf('b7','pt8','b7R',[1 5]); p.nc.lammax=5.1; p.sol.ds=0.1; p=cont(p,30); 
%% BD plot 
f=3; c=0; figure(f); clf; plotbra('b2R',f,c,'cl',p2pc('b1'),'lab',33); 
plotbra('b7R',f,c,'cl',p2pc('o2'),'lab',34); 
%% solns plot
torplot('b2R','pt33'); pause;  torplot('b7R','pt34'); 