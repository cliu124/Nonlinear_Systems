function [cj,aj,bj]=rbconv_Gjac(p,u) 
% jacobian for Rayleigh-Benard convection in vorticity formulation
par=u(p.nu+1:end); u=p.mat.fill*u(1:p.nu);
[ux,uz]=pdegrad(p.mesh.p,p.mesh.t,u); 
psix=ux(1,:); psiz=uz(1,:); 
omx =ux(2,:); omz =uz(2,:); 
Tx  =ux(3,:); Tz  =uz(3,:); 
R=par(1); 

cj=1; fu=zeros(p.nc.neq^2,1); fu(4)=-1; aj=-fu; 

bj=zeros(p.nc.neq^2*2,p.mesh.nt); 

bj(2*p.nc.neq*2+3,:)=R*ones(1,p.mesh.nt); % b231: (partial_x)_3 in second eqn.
bj(3,:)          =-omz;  % b211
bj(2*p.nc.neq+4,:)  =-psix; % b222 
bj(2*p.nc.neq+3,:)  =psiz;  % b221
bj(4,:)          =omx;   % b212

bj(5,:)          =ones(1,p.mesh.nt)-Tz;   % b311: (partial_x)_1 in third eqn.
bj(3*2*p.nc.neq,:)  =-psix; % b332
bj(3*2*p.nc.neq-1,:)=psiz;  % b331
bj(6,:)          =Tx;    % b312
