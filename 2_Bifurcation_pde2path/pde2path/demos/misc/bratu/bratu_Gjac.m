function [cj,afj,bj]=bratu_Gjac(p,u) 
% jacobian for Bratu 
par=u(p.nu+1:end); u=pdeintrp(p.mesh.p,p.mesh.t,u(1:p.nu)); 
fu=-10*(1-par(1)*exp(u)); cj=par(2); bj=0; afj=-fu; 