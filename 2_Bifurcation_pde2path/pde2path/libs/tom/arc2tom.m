function p=arc2tom(p) 
% arc2tom: convert arc-data-struc to tom (e.g., for meshref in t) 
%
%  p=arc2tom(p) 
% pad hopf.y with aux vars for pBC 
nu=p.nu; tl=p.hopf.tl; 
pad=p.hopf.y(1:nu,1)*ones(1,tl); 
p.hopf.y=[p.hopf.y; pad; p.hopf.T*ones(1,tl)]; 
par=p.u(p.nu+1:end); 
f=horhs(0,p.hopf.y(:,1),p,par); p.hopf.u0dot=f(1:p.nu); 
p.sw.para=3; p.sol.ds=p.sol.ds/5; % some safety for ds
p.hopf.lam=getlam(p)+p.sol.ds; %*p.hopf.tau(end); 
p.sol.restart=1; 
size(p.hopf.y), p.hopf.y(end,1)
return 