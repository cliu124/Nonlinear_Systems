function [c,a,f,b]=rbconv_G(p,u) 
% pde for Rayleigh-Benard convection in vorticity formulation
par=u(p.nu+1:end); u=p.mat.fill*u(1:p.nu);
[ux,uz]=pdegrad(p.mesh.p,p.mesh.t,u); 
u=pdeintrp(p.mesh.p,p.mesh.t,u);
om =u(2,:); 
% derivatives:
psix=ux(1,:); psiz=uz(1,:); 
omx =ux(2,:); omz =uz(2,:); 
Tx  =ux(3,:); Tz  =uz(3,:); 
R=par(1); c=1;
% variant 1: all in f
%b=0; a=0; 
% f1=-om; % Lap psi=om
%f2=R*P*Tx-(psix.*omz-psiz.*omx);
%f3=psix  -(psix.*Tz -psiz.*Tx);
% variant 2: use a, b
a=zeros(p.nc.neq^2,1); a(4)=1; % p.neq entries=diagonal a
b=zeros(p.nc.neq^2*2,1); % initialize; here 3^2*2=12
b(2*p.nc.neq*2+3)=R; % b231: (partial_x)_3 in second eqn.
b(2*2+1)=1;   % b311: (partial_x)_1 in third eqn.
f1=0*om;
f2=-(psix.*omz-psiz.*omx);
f3=-(psix.*Tz -psiz.*Tx);
f=[f1;f2;f3];
