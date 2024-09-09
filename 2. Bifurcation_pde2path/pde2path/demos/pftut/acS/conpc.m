function p=conpc(p,ds,n1,n2) % convenience function to continue with PC
% first n1 initial steps without PC, then n2 steps with PC 
p.sol.ds=ds; p=cont(p,n1); spplot(p); pause;
p.nc.nq=1; p.nc.ilam=[2 4]; p=cont(p,n2); spplot(p); 