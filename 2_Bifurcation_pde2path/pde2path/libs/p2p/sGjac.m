function J=sGjac(p,u)
% sGjac: dummy sGjac, set in stanparam; 
% useful to call jaccheck without ever setting sGjac 
fprintf('dummy sGjac from lib/p2p\n'); J=0*speye(p.nu); 