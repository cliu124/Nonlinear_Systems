function Gu=sGjacenn(p,u) % der.of PC missing 
Xbc=bcX(p,u); p.X(p.idx,:)=Xbc; % update BCs 
H0=u(p.nu+1); u=u(1:p.nu); N0=getN(p,p.X);
X=p.X+u.*N0; M=getM(p,X); LB=cotmatrix(X,p.tri); 
K=discrete_gaussian_curvature(X,p.tri); % same but faster
F=-1*spdiags(K,0,p.np,p.np)+2*H0^2*M; %H0 
Gu=-0.5*LB-F; Gu(p.idx,:)=0;
for i=1:length(p.idx); Gu(p.idx(i),p.idx(i))=1; end
