function p=conpc(p,ds,n1,n2) % convenience function to continue with PC
% first n1 initial steps without PC, then n2 steps with PC 
p.sol.ds=ds; p=cont(p,n1); %plotsol(p); pause;
p.fuha.qf=@qf2; p.fuha.qfder=@qf2der; 
p.nc.nq=2; p.nc.ilam=[6 7 5]; p=cont(p,n2);% plotsol(p); 