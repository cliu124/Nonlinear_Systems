function plotsolfu(dir,sfname,wnr,cnr,pstyle)
% adaption of plotsolf for SLOC 
p=loadp(dir,sfname); fprintf('lam=%g\n',getlam(p)); 
u=p.u; u(p.np+1:2*p.np)=-1./u(p.np+1:2*p.np); 
plotsolu(p,u,wnr,cnr,pstyle);  tname=[dir '/' sfname]; 
title(['v at ' tname]); end 
