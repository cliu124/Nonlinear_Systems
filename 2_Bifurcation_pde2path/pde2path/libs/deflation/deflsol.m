function [u,p]=deflsol(p,u1)  
% deflsol: solve deflated system 
p.fuha.sGb=p.fuha.sG; p.fuha.sGjacb=p.fuha.sGjac;  % mod rhs
p.fuha.sG=@deflsG; p.fuha.sGjac=@deflsGjac; 
%[Gu, Gun]=jaccheck(p);
[u,res,iter,Gu,Glam,p]=nloop(p,u1);
ga=deflfu(p,u);   
if res<p.nc.tol; % new soln found! 
   p.defl.nd=p.defl.nd+1; p.defl.u=[p.defl.u, u]; p.u=u; plotsol(p); 
   fprintf('sol found in deflation, res=%g, iter=%i, nd=%i, ga=%g\n', res,iter,p.defl.nd,ga); 
else fprintf('NO SOL in deflation, res=%g, iter=%i, ga=%g\n', res,iter,ga); 
end
p.fuha.sG=p.fuha.sGb; p.fuha.sGjac=p.fuha.sGjacb; 