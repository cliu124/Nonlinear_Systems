function qu=qjacV(p,u)
%par=u(p.nu+1:end); u=u(1:p.nu);%full domain
N0=getN(p,p.X); 
grX=grad(p.X,p.tri); grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:); grXz=grX(2*p.nt+1:3*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; grXz=E*grXz; %Iterpolated grad
%Translation
dphix=grXx*p.X; dphix=dot(dphix,N0,2); dphiy=grXy*p.X; dphiy=dot(dphiy,N0,2);
dphiz=grXz*p.X; dphiz=dot(dphiz,N0,2);
qux=(p.mat.fill'*dphix)'; quy=(p.mat.fill'*dphiy)'; quz=(p.mat.fill'*dphiz)';
qu=[qux; quy; quz]; M=getM(p,p.X); q4=sum(M,1); q4=(p.mat.drop*q4')'; 
size(qu), size(q4)
qu=[qu; q4];