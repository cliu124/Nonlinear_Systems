function bc=rbconv_bc_stressfree(p,u)
% stress free BC 
stifffac=10^4; % stiff spring factor
qd = stifffac*[[1;0;0] [0;1;0] [0;0;1]];% Dirichlet:
gv = zeros(p.nc.neq,1);% values
qm = qd; qm(3,3)=0; % mixed: Dirichlet for psi and om, and T_x=0

% ordering: bottom, right, top, left
bc = gnbc(p.nc.neq,qd,gv,qm,gv,qd,gv,qm,gv);