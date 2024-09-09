%% HP cont
p=hpcontini('s1','hpt1',2,'hpc1'); huclean(p); p.sw.bprint=2; 
%[Ja,Jn]=hpjaccheck(p); pause % to check the correct impl. of hpjac 
p.sol.ds=-0.01; p.nc.lammax=20; p.nc.lammin=-20; p.sw.jac=1; 
p.sw.spjac=0; p.nc.lammin=0.05; p.nc.tol=1e-3; p.nc.lammax=20; 
p.plot.bpcmp=4; p.nc.del=0.001; % seems needed 
%%
p=cont(p,10); 