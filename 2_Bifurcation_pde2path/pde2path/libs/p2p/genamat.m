function amat=genamat(p,Gu,Glam,tau,xi,xiq)
% GENAMAT: generate extended matrix for arclength system. 
%
%  amat=genamat(p,Gu,Glam,tau,xi,xiq)
% Gu, Glam = derivatives, xi=weighting of PDE variable vs. primary parameter, 
% xiq=weigthing of auxiliary equations.
%
% See also getder, xinorm
%xiq, xi, pause
if (p.nc.nq>0)
    amat=[[Gu Glam]; [xi*tau(1:p.nu)' xiq*tau(p.nu+1:p.nu+p.nc.nq)' (1-(xi+xiq))*tau(p.nu+p.nc.nq+1)]];  
else
    amat=[[Gu Glam]; [xi*tau(1:p.nu)' (1-xi).*tau(p.nu+p.nc.nq+1)]]; 
end
%figure(6); clf; spy(amat)