function [y1,r]=plotpp(p); 
u=p.cp.u; y1=u(3,:); 
r=sqrt(u(1,:).^2+u(2,:).^2); 
figure(10); plot(y1,r); 