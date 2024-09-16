function Gu=sGjac(p,u)  % PDE Jacobian for AC with pBC 
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
mu=par(1);
x1=up(1);
x2=up(2);

% Jacobian matrix of the nonlinear system
% \dot{x}_1=x1*[\mu-(x1^2+x2^2)]-x2
% \dot{x}_2=x2*[\mu-(x1^2+x2^2)]+x1
Gu(1,1)=mu-x1^2-x2^2-2*x1^2;
Gu(1,2)=-2*x2*x1-1;
Gu(2,1)=-2*x1*x2+1;
Gu(2,2)=mu-x1^2-x2^2-2*x2^2;

Gu=-Gu; %reverse the sign due to pde2path convention