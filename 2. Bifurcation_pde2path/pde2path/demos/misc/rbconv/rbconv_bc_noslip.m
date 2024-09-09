function bc=rbconv_bc_noslip(p,u)
% no slip BC 

stifffac=10^6; % stiff spring factor 
qd = stifffac*[[1;0;0] [0;1;0] [0;0;1]]; % Dirichlet 
qn = zeros(p.nc.neq); % zero flux:
gv = zeros(p.nc.neq,1); % values
% mixed: zero flux for psi and T, and psi=0
qm = qn; qm(2,1)=stifffac;
% ordering: bottom, right, top, left
bc = gnbc(p.nc.neq,qd,gv,qm,gv,qd,gv,qm,gv);