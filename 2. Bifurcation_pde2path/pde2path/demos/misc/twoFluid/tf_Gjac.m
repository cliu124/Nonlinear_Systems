function [cj,aj,bj]=tf_Gjac(p,u)  

% separate upde and auxiliary variables, here "par"
par=u(p.nu+1:length(u)); u=p.mat.fill*u(1:p.nu); 

% set parameters
delT=par(1); nu=par(2); L1=par(3); s=par(5);% - delT/L2*pi;

% linear part
cj=isoc([[nu,0,0];[0,nu,0];[0,0,1]],p.nc.neq,1); 

% one-vector, zero-vector
ov=ones(1,p.mesh.nt); zv=zeros(1,p.mesh.nt); 

bj=[
  % 1st eqn   2nd eqn  3rd eqn
    zv;(delT+s)*ov;   zv;zv;        zv;zv;   % 1st comp
    zv;zv;            zv;s*ov;      zv;zv;   % 2nd comp
    zv;-ov/L1;        zv;ov/L1;     zv;zv    % 3rd comp
];

aj=[0; 0; -1;
    0; 0; -1;
    0; 0;  0];

% nonlinear part

% gradients
[ux1,ux2]=pdegrad(p.mesh.p,p.mesh.t,u); 
ux1=ux1; ux2=ux2;
u1x1=ux1(1,:); u1x2=ux2(1,:);
u2x1=ux1(2,:); u2x2=ux2(2,:);
u3x1=ux1(3,:); u3x2=ux2(3,:);

bj = bj + [
     u3x2;-u3x1;    zv;zv;      zv;zv;
       zv;zv;     u3x2;-u3x1;   zv;zv;
    -u1x2;u1x1;  -u2x2;u2x1;    zv;zv
];

end