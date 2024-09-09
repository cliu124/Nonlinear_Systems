para=4; ds=0.1; figure(2); clf;  aux=[]; %aux.tl=60; 
for m=[20,40,60]
  aux.tl=m; 
  p=hoswibra('hom1d','hpt1',ds,para,['m' mat2str(m)],aux); nsteps=10; %20+m; 
  p.hopf.jac=1; p.nc.dsmax=0.5; p.hopf.xi=0.05; p.file.smod=5; p.sw.verb=0; 
  p.hopf.flcheck=1; % switch for Floquet-comp: 0: off, 1:floq, 2: floqps 
  bw=1; beltol=1e-6; belimax=5; % border-width, bel-parameters 
  p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
  t1=tic; p=cont(p,nsteps); toc(t1) 
end 
%% BD T for convergence illustration 
bpcmp=6; wnr=3; figure(wnr); clf; 
plotbra('1db1','pt30',wnr,bpcmp,'lsw',0); 
plotbra('m40','pt60',wnr,bpcmp,'cl','r','lsw',0,'tyst','--'); 
plotbra('m60','pt80',wnr,bpcmp,'cl','b','lsw',0,'tyst','-.'); 
axis([0.5 1 7.3 7.5]); 
xlabel('r'); ylabel('T'); 
% add comparison to analytical soln 
p=loadp('hom1d','hpt1'); k2=0; rstep=0.075; ms=10; plotana1(k2,rstep,'k*',2,ms); 
box on; 