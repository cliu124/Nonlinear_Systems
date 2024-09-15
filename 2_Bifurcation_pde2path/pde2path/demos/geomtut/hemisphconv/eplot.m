function ev=eplot(p)
n=p.np; %par=p.u(p.nu+1:end); 
r1=pderesi(p,p.u); mr1=max(abs(r1));  %p.up(1:n)=res; pplot(p,10); title('r'); 
X=p.X; T=p.tri; Xn=normrow(X); Xd=abs(Xn-1); mXd=max(Xd); 
p.sw.msw=0; M=getM(p,X); p.sw.msw=1; Mf=getM(p,X); 
l2=sum(M*(Xd.^2),1); l2=sqrt(l2); 
LB=cotmatrix(X,T); N=getN(p,X); H2=0.5*dot(LB*X,N,2); H=M\H2; Hf=Mf\H2;  
[k,H1,K1]=discrete_curvatures(X,T); K=M\K1; Kf=Mf\K1; 
id=p.idx; bulk=setdiff(1:n,id); 
Hd=H(bulk)+1; Hfd=Hf(bulk)+1; 
Kd=K(bulk)-1; Kfd=Kf(bulk)-1; 
Hdm=max(abs(Hd)); Hfdm=max(abs(Hfd)); 
Kdm=max(abs(Kd)); Kfdm=max(abs(Kfd)); 
Mb=M(bulk,bulk); Mfb=Mf(bulk,bulk);
Hdl2=sqrt(sum(Mb*(Hd.^2))); Hfdl2=sqrt(sum(Mfb*(Hfd.^2)));
Kdl2=sqrt(sum(Mb*(Kd.^2))); Kfdl2=sqrt(sum(Mfb*(Kfd.^2)));
ev=[n; mr1; Hdm; Hfdm; Kdm; Kfdm; Hdl2; Hfdl2; Kdl2; Kfdl2]; 
%   1   2   3    4      5    6     7     8      9     10 
% solve first with M_Vor 
p.sw.msw=0; [p,res]=solfixpar(p); res
X=p.X; Xn=normrow(X); Xd=abs(Xn-1); mXd=max(Xd); 
M=getM(p,X); 
l2=sum(M*(Xd.^2),1); l2=sqrt(l2); 
ev=[ev; res; mXd; l2]; 
%        11   12   13 
p.sw.msw=1; [p,res]=solfixpar(p); res
X=p.X; Xn=normrow(X); Xd=abs(Xn-1); mXd=max(Xd); 
M=getM(p,X); 
l2=sum(M*(Xd.^2),1); l2=sqrt(l2); 
ev=[ev; res; mXd; l2]; 
%        14   15  16
A=doublearea(X,T)/2; e=edge_lengths(X,T);
el=max(e,[],2); % long edges 
s=sum(e,2)/2; r=A./s; q=max(el./r); % inradii.and "quality" 
% R=e(:,1).*e(:,2).*e(:,3)./(8*A);  % outradius
a1=max(A); a2=min(A); 
h1=max(el); h2=min(el); 
vl=vallist(p); val=mean(vl); 
ev=[ev; q;   a1; a2; h1; h2; val]; 
%       17   18      20




