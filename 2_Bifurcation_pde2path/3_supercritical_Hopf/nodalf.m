function f=nodalf(p,u) % nonlinearity for cGL 
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
mu=par(1);
x1=up(1);
x2=up(2);
f(1,1)=x1*(mu-x1^2-x2^2)-x2;
f(2,1)=x2*(mu-x1^2-x2^2)+x1;

f=-f;
