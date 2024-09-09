function r=lvsG(p,u)  % rhs for lvoc. 
n=p.np; par=u(p.nu+1:end); beta=par(1); d1=par(4); d2=par(5); 
ga1=par(6); ga2=par(7); rho=par(8); p1=par(11); p2=par(12); % extract pars 
v1=u(1:n); v2=u(n+1:2*n); l1=u(2*n+1:3*n); l2=u(3*n+1:4*n); % extract fields 

f1=v1.*(1-beta*v1-v2); f2=(v1-1).*v2; % bulk nonlin. 
a11=1-2*beta*v1-v2; a12=-v1; a21=v2; a22=v1-1; 
f3=rho*l1-a11.*l1-a21.*l2; 
f4=rho*l2-a12.*l1-a22.*l2; 
f=[f1; f2; f3; f4];  F=p.mat.M*f; 

l11=l1(1); l21=l2(1); [h,hv,hk]=hfu(p,u); % bd vals of lambda, and harvest 
g1=-ga1*h(1); g2=-ga2*h(2); % BC for states 
g3=-(p1-ga1*l11)*hv(1); g4=-(p2-ga2*l21)*hv(4); % BC for co-states 

F(1)=F(1)+g1; F(n+1)=F(n+1)+g2; % add BC directly into F, states 
F(2*n+1)=F(2*n+1)+g3; F(3*n+1)=F(3*n+1)+g4; % BC for co-states 
zM=0*speye(n); K=p.mat.K; % finally assemble the system K 
K=[[d1*K zM zM zM]; [zM d2*K zM zM]; [zM zM -d1*K zM]; [zM zM zM -d2*K]]; 
r=K*[v1;v2;l1;l2]-F;  % and compute the residual as usual 