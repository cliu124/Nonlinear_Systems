function [L,M,N]=get2ff(p,X) % 2nd fundamental form 
Xxx=p.mat.Kx*(p.mat.Dx*X); Xxy=p.mat.Kx*(p.mat.Dy*X);Xyy=p.mat.Ky*(p.mat.Dy*X);
nu=getN(p,X);  L=dot(nu,Xxx,2); M=dot(nu,Xxy,2); N=dot(nu,Xyy,2); 