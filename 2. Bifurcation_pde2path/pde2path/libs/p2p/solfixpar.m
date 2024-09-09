function [p,res,Gu]=solfixpar(p)
% solhere: run Newtonloop to solve at fixed pars
[u1,res,iter,Gu,Glam,p]=nloop(p,p.u); 
if res<p.nc.tol 
    if p.sw.Xcont>0; [p,u1]=updX(p,u1);  end
    p.u=u1;
end
%r1=pderesi(p,p.u); mr2=max(abs(r1)); mr2