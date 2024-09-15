function [jcaval,jcav,jcavd]=jcaiT(p,sol,rho) 
% jcai: jca along canonical path, with dependence
T=sol.par(1); tv=T*sol.t; tl=length(tv); % T, tv(end), pause 
jcav=zeros(1,tl); par=p.u(p.nu+p.nc.nq+1:end); 
for i=1:tl
    jcav(i)=jca(p,[sol.u(1:p.nu+p.nc.nq,i);par]);
    jcavd(i)=exp(-rho*tv(i))*jcav(i); %jca(p,[sol.u(1:p.nu+p.nc.nq,i);par]);
end
jcaval=trapz(tv,jcavd); 
