function Gu=sGnodDjac(p,u) % der.of PC missing 
u=u(1:p.nu); N0=getN(p,p.X);
X=p.X+u.*N0; M=getM(p,X); LB=cotmatrix(X,p.tri); 
[~,Hr,Kr]=discrete_curvatures(X,p.tri); Hr=M\Hr;Kr=M\Kr;  F=4*Hr.^2-2*Kr;
F=spdiags(F,0,p.np,p.np);
Gu=-LB-M*F; Gu(p.idx,:)=0;
for i=1:length(p.idx); Gu(p.idx(i),p.idx(i))=p.sf; end
