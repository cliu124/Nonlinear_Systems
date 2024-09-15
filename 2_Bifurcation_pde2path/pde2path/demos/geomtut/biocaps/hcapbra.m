function out=hcapbra(p,u) % mod of cmcbra to put Helfrich energy E on branch 
par=u(p.nu+1:end); al=par(1); l1=par(2); c0=par(3); b=par(4);  
id=p.idx; X=p.X; V=getV(p,u); A=getA(p,u); 
M=getM(p); M=M(1:p.np,1:p.np); H=p.u(p.np+1:p.nu); kapg=zeros(p.np,1); 
N1=getN(p,X); c=cross(N1(id,:),-X(id,:)/al,2); 
kapg(id)=-sign(N1(id,3)).*vecnorm(c,2,2); % geodesic curvature 
bdE1=2*pi-bdint(p,kapg); % by kap_g
K=discrete_gaussian_curvature(X,p.tri); K(id)=0; bdE2=sum(K); % for comparison 
E=sum(M*(H-c0).^2); 
E1=E+l1*A+b*bdE1;
E2=E+l1*A+b*bdE2; 
mqd=meshqdat(p);
out=[par; V; A; E1; E2; bdE1;mqd]; 
%   1-5   6  7  8   9   10   11