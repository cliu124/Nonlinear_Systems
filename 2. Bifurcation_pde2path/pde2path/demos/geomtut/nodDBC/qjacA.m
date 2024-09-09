function qu=qjacA(p,u) % changed cause sign of N flipped 
par=u(p.np+1:end); u=u(1:p.np); N=getN(p,p.X); 
X=p.X+u.*N; M=getM(p,X); qu=-2*par(1)*sum(M);