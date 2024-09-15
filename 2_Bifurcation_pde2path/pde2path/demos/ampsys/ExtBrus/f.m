function w=f(u,p)
a=p.par(1);b=p.par(2);
f1=a-(b+1)*u(1)+u(1)^2*u(2)-u(1)+u(3);
f2=b*u(1)-u(1)^2*u(2);
f3=u(1)-u(3);
w=[f1;f2;f3];