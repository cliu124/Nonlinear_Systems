function p=oosetfemops(p) 
M=massmatrix(p.X,p.tri,'full'); p.mat.M=filltrafo(p,M); 
tol=1e-3; X=p.X;
m1=min(X(:,1)); M1=max(X(:,1)); p.idxx=find(abs(X(:,1)-m1)<tol |abs(X(:,1)-M1)<tol);
m2=min(X(:,2)); M2=max(X(:,2)); p.idxy=find(abs(X(:,2)-m2)<tol | abs(X(:,2)-M2)<tol);
m3=min(X(:,3)); M3=max(X(:,3)); p.idxz=find(abs(X(:,3)-M3)<tol | abs(X(:,3)-m3)<tol);
p.idxz1=find(abs(X(:,3)-M3)<tol ); p.idxz2=find(abs(X(:,3)-m3)<tol );
p.idxx=setdiff(p.idxx,p.idxz);p.idxy=setdiff(p.idxy,p.idxz);
X=p.mat.drop*p.X; 
%m3=min(X(:,3)); M3=max(X(:,3)); p.idxz1=find(abs(X(:,3)-M3)<tol);
 %p.idxz2=find(abs(X(:,3)-m3)<tol);
p.idxzd=find(abs(X(:,3)-M3)<tol | abs(X(:,3)-m3)<tol);
end


