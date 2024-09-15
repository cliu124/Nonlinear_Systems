function r=sGbdcurve(p,u)  % PDE rhs, discrete mean curvature 
if p.bcsw==2; Xbc=bcX(p,u); p.X(p.idx,:)=Xbc; end % set Bdcurve 
par=u(p.nu+1:end); H0=par(1); N=getN(p,p.X); X=p.X+u(1:p.np).*N; 
M=getM(p,X); LB=cotmatrix(X,p.tri); N=getN(p,X); 
r=-0.5*dot(LB*X,N,2)+M*(H0*ones(p.np,1)); % PDE-rhs, i.e., H-H0=0
switch p.bcsw   % BCs  
    case 1; % X_3=al*sin(k*phi) 
        al=par(4); k=par(5); phi=angle(X(p.idx,1)+1i*X(p.idx,2)); 
        r(p.idx)=X(p.idx,3)-al*sin(k*phi);
    otherwise; r(p.idx)=u(p.idx);  % \pa X=\ga  (boundary curve) 
end 