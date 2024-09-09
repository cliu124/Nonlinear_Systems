close all; format compact; keep pphome; 
%% Bif to standing waves, sometimes must be enforced by phase cond. 
para=4; ds=0.2; dsmax=1.1; xi=1e-2; nsteps=10; 
figure(2); clf; aux=[]; aux.tl=15; aux.z=[1 1i]; tol=1e-6; 
pc=1; % run with (pc=1) or without (pc=0) rotational PC. pc=1 forces SW! 
if pc; aux.nqh=1; aux.qfh=@qfh; aux.qfhder=@qfhjac; end 
for bnr=6; %[2 3 5 6]
switch bnr
 case 2; p=hoswibra('tr','hpt2',ds,para,'sw2',aux); 
 case 3; p=hoswibra('tr','hpt3',ds,para,'sw3',aux); 
 case 5; aux.z=[1 -4i]; p=hoswibra('tr','hpt5',ds,para,'sw5',aux); 
 case 6; aux.z=[1 8i]; tol=1e-6; dsmax=0.5; p=hoswibra('tr','hpt6',ds,para,'sw6',aux); 
end 
p.hopf.xi=xi; p.hopf.jac=1; p.nc.dsmax=dsmax; p.usrlam=[]; p.nc.tol=tol; 
p.u(p.nu+2:p.nu+3)=[0;1]; % default case 
p.hopf.flcheck=0; p.file.smod=1; % switch off floquet (if somewhat slow)
bw=1; if pc; bw=2; end 
AMG=0; p.sw.verb=3;  % set AMG=1 if ilupack available
if ~AMG; p=setbel(p,bw,1e-3,10,@lss); % use BEL without ilupack 
else p=setilup(p,1e-3,200); p=setbel(p,bw,1e-3,10,@lssAMG); end
if pc; p.nc.ilam=1; p.hopf.ilam=4; p.u0x=5*p.mat.Krot*p.hopf.tau(1:p.nu)'; end 
pause
tic; p=cont(p,nsteps);  toc
end 
%% BD, L^2
bpcmp=5; figure(3); clf; plotbra('tr',3,bpcmp,'cl','k'); 
plotbra('sw2','pt20',3,bpcmp,'cl','r'); plotbra('sw3','pt10',3,bpcmp,'cl','m'); 
plotbra('sw4','pt10',3,bpcmp,'cl','b'); plotbra('sw5','pt10',3,bpcmp,'cl','r'); 
plotbra('sw6','pt10',3,bpcmp,'cl','m'); 
plotbra('h2','pt10',3,bpcmp,'cl','r','tyun','--'); 
plotbra('h3','pt10',3,bpcmp,'cl','m','tyun','--'); 
plotbra('h5','pt10',3,bpcmp,'cl','b','tyun','--'); 
plotbra('h6','pt10',3,bpcmp,'cl','r','tyun','--'); 
plotbra('h7','pt10',3,bpcmp,'cl','m','tyun','--'); 
axis([-0.2 0.6 0 0.9]); xlabel('r'); ylabel('||u||_*'); box on; 
%% soln plot, here use specialized 'hoplot' ! 
pstyle=2;
hoplotrotf('sw2','pt20',4,1,pstyle); pause; 
hoplotrotf('sw3','pt4',4,1,pstyle); pause; 
hoplotrotf('sw4','pt10',4,1,pstyle); pause; 
hoplotrotf('sw5','pt10',4,1,pstyle);   
%% Floquet a posteriori, choose points of interest! 
aux.nfloq=20; [muv1,~,~,h]=floqpsap('sw3','pt5',aux); %[muv1,~,~,~,h]=floqap('sw3','pt5',aux); 