%% HP and BP continuations, first the wave (Turing-Hopf) branch  
p=hpcontini('hom1d','hpt1',1,'hpc1'); huclean(p); p.sw.bprint=2;p.plot.bpcmp=2; 
%[Ja,Jn]=hpjaccheck(p); pause % to check the correct impl. of hpjac 
p.sol.ds=-0.01; p.nc.lammax=1.25; p.nc.lammin=0.8; p0=p; p=cont(p,20); 
p=p0; p=setfn(p,'hpc1b'); p.sol.ds=-p.sol.ds; p=cont(p,20); % other direction 
%% check HP continuation 
p=hpcontexit('hpc1','pt5','t1'); % puts the HBP into dir 't1' 
p=hoswibra('t1','hpt1',ds,para,'t1h'); p.nc.lammax=3.5; p=cont(p,5); % continue
%% BP-cont for Turing: 
p=bpcontini('hom1d','bpt1',1,'bpc1'); p.sw.bprint=2; p.plot.bpcmp=2; 
p.sol.ds=-0.01; p.nc.lammin=0.8; p.nc.lammax=1.15; p0=p; p=cont(p,20); 
p=p0; p=setfn(p,'bpc1b'); p.sol.ds=-p.sol.ds; p=cont(p,20); % other direction 
%% check BP continuation 
p=bpcontexit('bpc1','pt5','t2'); % puts the BP into dir 't2'
p=swibra('t2','bpt1','t2b'); p.nc.lammax=3.5; p=cont(p,5); % swibra and cont
%% HP-cont for 'hopf'.branch 
p=hpcontini('hom1d','hpt2',1,'hpc2'); p.sw.bprint=2; p.plot.bpcmp=2; 
p.sol.ds=-0.01; p.nc.lammin=0.8; p.nc.lammax=1.25; p0=p; p=cont(p,20); 
p=p0; p=setfn(p,'hpc2b'); p.sol.ds=-p.sol.ds; p=cont(p,20); % other direction 
%% plot BD of b over a 
bpcmp=2; figure(3); clf; 
plotbra('hpc1',3,bpcmp,'cl','r','lab',3,'lms',0.1); plotbra('hpc1b',3,bpcmp,'cl','r'); 
plotbra('hpc2',3,bpcmp,'cl','m'); plotbra('hpc2b',3,bpcmp,'cl','m','lab',5,'lms',0.1); 
plotbra('bpc1',3,bpcmp,'cl','b','lab',3,'lms',0.1); plotbra('bpc1b',3,bpcmp,'cl','b'); 
axis([0.8 1.1 2.4 3.3]); xlabel('a'); ylabel('b'); 