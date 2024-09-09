function r=resi(p,u)
% RESI: residual for PDE and auxiliary equations
%
%  r=resi(p,u)
%
% See also pderesi, getGu
switch p.sw.spcont
    case 0; r=pderesi(p,u);  % regular continuation, pde part              
      if(p.nc.nq>0) r=[r;p.fuha.qf(p,u)]; end  % possibly add q-part 
    case 1; r=bpcresi(p,u); % BP continuation       
    case 2; r=fpcresi(p,u); % fold point continuation 
    case 3; r=hpcresi(p,u); % HP continuation     
end