function qu=qjacV(p,u)
% qjacV: u-derivative of volume constraint qfV
M=getM(p,p.X); qu=sum(M,1);
%u=u(1:p.np); N=getN(p,p.X); X=p.X+u.*N; M=getM(p,X); qu=sum(M,1);