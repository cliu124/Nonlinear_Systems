function ja=jca(p,u) % J_c average, overload here since J_c is on boundary 
% (average makes no sense) 
ja=p.fuha.jcf(p,u); 