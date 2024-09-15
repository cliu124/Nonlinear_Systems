function qu=qfjac(p,u)  
M=p.mat.M(1:p.np,1:p.np); 
qu=([M*ones(p.np,1); M*ones(p.np,1)])'/p.vol; 