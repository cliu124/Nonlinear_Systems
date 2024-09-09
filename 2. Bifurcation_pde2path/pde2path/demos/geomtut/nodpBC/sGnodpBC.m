function r=sGnodpBC(p,u)  % nodoids with pBCs 
par=u(p.nu+1:end); H0=par(1); del=par(6); % height 
sx=par(7); sy=par(8); sz=par(9); sphi=par(10);% Lag.for translations & rotation 
u=u(1:p.nu); uf=p.mat.fill*u; N0=getN(p,p.X); X=p.X+uf.*N0; % fill, then as usual
m=min(X(p.idx,3)); M=max(X(p.idx,3)); old=M-m; % old height 
scaz=del/old;  a=[ones(p.np,2) scaz*ones(p.np,1)]; X=a.*X; % scale X in z-direction 
M0=getMf(p,X); M=p.mat.fill'*M0; % like fill-trafo (fm'*M0*fm) 
LB=cotmatrix(X,p.tri); LB=p.mat.fill'*LB; N=getN(p,X); 
%Derivatives of phaseconditions at p.X
grX=grad(p.X,p.tri); grXx=grX(1:p.nt,:); 
grXy=grX(p.nt+1:2*p.nt,:); grXz=grX(2*p.nt+1:3*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; grXz=E*grXz;
dphix=grXx*p.X; dphix=dot(dphix,N0,2); 
dphiy=grXy*p.X; dphiy=dot(dphiy,N0,2);
dphiz=grXz*p.X; dphiz=dot(dphiz,N0,2);
dphi=(-p.X(:,2).*grXx*p.X+p.X(:,1).*grXy*p.X)./(p.X(:,1).^2+p.X(:,2).^2); 
dphi=dot(dphi,N0,2);
r=dot(-LB*X-M*(2*H0*N),p.mat.fill'*N,2)+M*(sx*dphix+sy*dphiy+sz*dphiz+sphi*dphi); 