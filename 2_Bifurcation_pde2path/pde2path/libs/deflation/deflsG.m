function r=deflsG(p,u)  
% deflsG: sG modified for deflation 
% see deflinit for meaning of parameters
global p2pglob; 
r1=p.fuha.sGb(p,u); ga=deflfu(p,u); r=ga*r1; 
p2pglob.cvec=r1; 