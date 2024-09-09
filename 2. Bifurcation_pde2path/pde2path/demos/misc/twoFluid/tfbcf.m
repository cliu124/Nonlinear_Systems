function [c,a,f,b]=tff(p,u)  

% separate upde and auxiliary variables, here "par"
par=u(p.nu+1:end); u=p.fill*u(1:p.nu); 

% set parameters
delT=par(1); nu=par(2); L1=par(3); s=par(5);% + delT/L1*pi; %deviation from onset speed

% pde
c=isoc([[nu,0,0];[0,nu,0];[0,0,1]],p.asw.neq,1); 

b=[
  % 1st eqn   2nd eqn  3rd eqn
    0;delT+s; 0;0;     0;0;   % 1st comp
    0;0;      0;s;     0;0;   % 2nd comp
    0;-1/L1;   0;1/L1; 0;0  % 3rd comp
];

a=[
    0; 0; -1;
    0; 0; -1;
    0; 0;  0
];

% gradients
[ux1,ux2]=pdegrad(p.mesh.p,p.mesh.t,u); 
u1x1=ux1(1,:); u1x2=ux2(1,:);
u2x1=ux1(2,:); u2x2=ux2(2,:);
E1=-ux1(3,:); E2=-ux2(3,:);

f1=( - (E2.*u1x1-E1.*u1x2) );
f2=( - (E2.*u2x1-E1.*u2x2) );
f3=zeros(1,p.mesh.nt);

f=[f1;f2;f3];
%f=[0;0;0];