function qu=fchqjac(p,u)
M=p.mat.Ms; n=p.nu/2; 
qu=([M*ones(n,1)/p.vol; 0*ones(n,1)])';
%qu=([M*ones(n,1); 0*ones(n,1)])';


