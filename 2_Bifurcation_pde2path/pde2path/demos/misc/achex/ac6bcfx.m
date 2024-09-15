function bc=ac6bcfx(p,u) 
% BC for geometry with six edges, mixed Neumann and Dirichlet, x-dependent 
sf=10^4; % stiff spring factor for Dirichlet part
qd = mat2str(sf);gd = [mat2str(sf*p.u(p.np+1)) '*x']; % pseudo-Dirichlet: 
qn = '0'; gn = '0'; % Neumann part
% ordering in geometry from rec: bottom, right, top, left
bc = gnbcs(p.nc.neq,qn,gn,qn,gn,qn,gn,qn,gn,qd,gd,qn,gn);


