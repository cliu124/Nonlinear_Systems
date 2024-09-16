function r=sG(p,u)  % AC with periodic BC 
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
mu=par(1);
x1=up(1);
x2=up(2);

% Right hand side of 
% \dot{x}_1=x2
% \dot{x}_2=\mu*(x1+x2)-x2-x1^3-3*x1^2*x2
r(1,1)=x2;
r(2,1)=mu*(x1+x2)-x2-x1^3-3*x1^2*x2;

r=-r; %reverse the sign due to pde2path convention