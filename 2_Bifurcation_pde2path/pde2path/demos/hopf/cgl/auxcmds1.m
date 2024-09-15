%% C1: change continuation param 
p=hoswiparf('1db1','pt39','c5b',5,0.1); clf(2); p.usrlam=0.25; p=cont(p,20);
%% C2: plot BD and solns 
bpcmp=6; pstyle=3; wnr=3; figure(wnr); clf; 
plotbra('c5b','pt19',3,bpcmp,'lsw',0); xlabel('c_5'); ylabel('T'); 
p=loadp('c5b','pt20'); plotana2(p,0,5,0.25,1,0.05,'b*',2,1); axis tight; 
hoplotf('1db1','pt28',1,1); figure(1); title('c_5=1'); 
%% C3: mesh-refinement in t using TOM:  
p=loadp('1db2','pt10','1db1ref'); hoplot(p,4,1,1); 
% switch to nat.-parametr., and reset tolerances, then cont 
p=arc2tom(p); p.hopf.tau=[]; p.sol.ds=0.01; 
p.hopf.tom.AbsTol=5e-5; p.hopf.tom.RelTol=5e-4; p=cont(p,5); 
p.sw.para=4; p=tom2arc(p); p=cont(p,5); % switch back to arclength and cont
%% C4: mesh-refinement using hopftref
p=loadp('1db3','pt17','1db1ref'); 
% hogradinf(p); % info about max_t |udot| (here useless, since u is harmonic) 
p=hopftref(p,4); p=cont(p,5); % bisect 5 intervals after t=4 and cont 
%% C5: uniform mesh-refinement 
p=loadp('1db3','pt17','1db1ref'); fac=2.3; 
p=uhopftref(p,fac); p=cont(p,5); % increase # time-slices by fac, then cont  