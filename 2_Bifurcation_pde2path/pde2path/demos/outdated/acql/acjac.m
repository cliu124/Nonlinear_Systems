function [cj,aj,bj]=acjac(p,u)  % jacobian for AC quasilin 
% -div[(al+del*u+ga*u^2)*grad u]-lam u-u^3+u^5=0
par=u(p.nu+1:end); al=par(1); del=par(2); ga=par(3); lam=par(4); 
u=u(1:p.nu); [po,t,e]=getpte(p); nt=size(t,2); 
[ux,uy]=pdegrad(po,t,u); 
uxn=pdeprtni(po,t,ux); % tria2nodal interpol 
uyn=pdeprtni(po,t,uy);
[uxx,uxy]=pdegrad(po,t,uxn); 
[uxy,uyy]=pdegrad(po,t,uyn); 
lapu=uxx+uyy; % defined on triangles 
ugradsq=ux.*ux+uy.*uy; 
u=pdeintrp(po,t,u); 
cj=al+del*u+ga*u.^2; 
aj=-(lam+3*u.^2-5*u.^4+del*lapu+2*ga*(u.*lapu+ugradsq)); 
bj=zeros(p.nc.neq*p.nc.neq*2,nt); 
bj(1,:)=del*ux+2*ga*u.*ux; 
bj(2,:)=del*uy+2*ga*u.*uy;
