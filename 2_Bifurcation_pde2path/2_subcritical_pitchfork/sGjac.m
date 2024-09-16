function Gu=sGjac(p,u)  % PDE Jacobian for AC with pBC 
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
mu=par(1);
x1=up(1);
x2=up(2);
Gu(1,1)=mu+3*x1^2-5*x1^4;
Gu(1,2)=0;
Gu(2,1)=0;
Gu(2,2)=-1;

Gu=-Gu; %reverse the sign due to pde2path convention