function r=resiGuphi(t,au)
% resiGuphi: interface for numjac for \pa_u(Gu*phi) 
% currently only works with p.sw.jac=1; 
global pGu % the parameters and psi from pGu
p=pGu;
ua=[au; p.u(p.nu+1:end)]; phi=p.phi; % p.u(p.nu+1:end)', pause 
r0=pderesi(p,ua); Gu=getGupde(p,ua,r0); r=Gu*phi; 
%p.sw.jac=1; Gu2=getGupde(p,ua,r0); Gd=abs(Gu-Gu2); 7, full(Gd(1:5,1:5)), 
%max(max(Gd)), pause 

