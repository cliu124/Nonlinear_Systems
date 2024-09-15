function N=getN(p,X)  
n=cross(p.mat.Dx*X,p.mat.Dy*X,2); nn=sqrt(n(:,1).^2+n(:,2).^2+n(:,3).^2); N=n./nn;    