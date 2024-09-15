function y=v2tom(p,yv) 
% v2tom: convert soln vector to matrix for tom
nu=p.nu+p.nc.nq; y=reshape(yv(1:nu*p.hopf.tl),nu,p.hopf.tl); 