function J=mynumjac(f,v,r); % for testing numjac
global pGu; 
p=pGu;  %p.phi', pause 
nu=p.nu; del=p.nc.del; J=sparse(nu,nu); 
for j=1:nu
  up=v+p.nc.del*ej(j,nu); 
  r1=feval(f,0,up); 
  r2=(r1-r)/del;  
  %J(:,j)=J(:,j)+r3; 
  J=J+sparse(1:nu,j,r2,nu,nu); 
end