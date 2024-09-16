function Gu=sGjac(p,u)  % PDE Jacobian for AC with pBC 
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
mu=par(1);
x1=up(1);
x2=up(2);

% Gu is Jacobian matrix associated with nonlinear system:
% \dot{x}_1=x2
% \dot{x}_2=\mu*(x1+x2)-x2-x1^3-3*x1^2*x2
Gu(1,1)=0;
Gu(1,2)=1;
Gu(2,1)=mu-3*x1^2-6*x1*x2;
Gu(2,2)=mu-1-3*x1^2;

Gu=-Gu; %reverse the sign due to pde2path convention