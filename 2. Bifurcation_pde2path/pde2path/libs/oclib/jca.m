function ja=jca(p,u) 
% jca: integrate J_c in space and normalize by p.vol 
jc=p.fuha.jcf(p,u); 
ja=sum(p.mat.M(1:p.np,1:p.np)*jc)/p.vol; 