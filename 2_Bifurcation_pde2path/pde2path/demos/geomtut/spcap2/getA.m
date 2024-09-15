function A=getA(p,u)
u=u(1:p.nu); N0=getN(p,p.X); X=p.X+u.*N0;
n=cross(p.mat.Dx*X,p.mat.Dy*X,2);
A=sum(p.mat.M*(sqrt(dot(n,n,2))));