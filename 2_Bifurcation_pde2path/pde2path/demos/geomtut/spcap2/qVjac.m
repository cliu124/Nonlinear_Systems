function qu=qVjac(p,u) % (approximate) derivative of volume 
N0=getN(p,p.X); X=p.X+u(1:p.nu).*N0; 
N=cross(p.mat.Dx*X,p.mat.Dy*X,2); % normal at X, NOT normalized 
N2=dot(N0,N,2); % volume element (in direction u) since M lives in preimage 
MX=p.mat.M*spdiags(N2,0,p.nu,p.nu); 
qu=sum(MX,1); 