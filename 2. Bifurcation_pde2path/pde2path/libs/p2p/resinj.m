function r=resinj(t,au)
% RESINJ: interface for numjac for pderesi
global pj % get the parameters from pj 
p=pj; ua=[au;p.u(p.nu+1:end)];  r=pderesi(p,ua); 
%12, ua(1:4)', ua(p.nu+1:end)',r(1:4)', pause 