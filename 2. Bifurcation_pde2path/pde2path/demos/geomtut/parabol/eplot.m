function ev=eplot(p)
n=p.np; par=p.u(p.nu+1:end); a=par(1); b=par(2); 
r1=pderesi(p,p.u); mr1=max(r1); % resi 
X=p.X; T=p.tri; x=X(:,1); y=X(:,2); z=X(:,3); 
za=zfu(p,x,y); r2=z-za; mr2=max(abs(r2)); % zfu=exact z
p.up(1:n)=r2; pplot(p,10); title('z-z_a'); nola
mrbd=max(abs(r2(p.idx))); 
M=getM(p,X); l2=sum(M*(r2.^2),1); l2=sqrt(l2); 
LB=cotmatrix(X,T); N=getN(p,X); 
a2=a^2; b2=b^2; a4=a2^2; b4=b2^2; 
Hf=(a2+b2+4*x.^2/a2+4*y.^2/b2)./(a2*b2*(1+4*x.^2/a4+4*y.^2/b4).^1.5); % exact H 
Hn=-0.5*dot(LB*X,N,2); 
id=p.idx; bu=setdiff(1:n,id); Mb=M(bu,bu); % bulk nodes 
Hnn=M\Hn; Hd1=Hnn(bu)+Hf(bu); mh1=max(abs(Hd1)); 
hl21=sum(Mb*(Hd1.^2),1); hl21=sqrt(hl21); 
Hd2=Hn(bu)-Mb*Hf(bu); mh2=max(abs(Hd2));
l=edge_lengths(X,T); l=l(:); maxl=max(l); minl=min(l); 
ev=[p.np; maxl; minl; mr1;  mr2;   l2;  mh1;  hl21; mh2; mrbd;]; 
%    1                 4     5      6    7     8     9
%                     resi  z-err       H-err      weak-H
A=doublearea(X,T)/2; e=edge_lengths(X,T);
el=max(e,[],2); % long edges 
s=sum(e,2)/2; r=A./s; q=max(el./r); % inradii.and "quality" 
% R=e(:,1).*e(:,2).*e(:,3)./(8*A);  % outradius
a1=max(A); a2=min(A); 
h1=max(el); h2=min(el); 
vl=vallist(p); val=mean(vl); 
ev=[ev; q;   a1; a2; h1; h2; val]; 
%       11                   16 