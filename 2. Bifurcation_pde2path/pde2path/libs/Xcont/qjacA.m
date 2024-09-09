function qu=qjacA(p,u)
% qjacA: u-derivative of area constraint qfA
par=u(p.np+1:end); u=u(1:p.np); N=getN(p,p.X); 
X=p.X+u.*N; M=getM(p,X); qu=2*par(1)*sum(M);