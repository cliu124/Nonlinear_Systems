function A=getA(p,u) 
u=u(1:p.nu); N=getN(p,p.X); X=p.X+(p.mat.fill*u).*N; 
Af=doublearea(X,p.tri); A=sum(Af); %  check sum(M) 