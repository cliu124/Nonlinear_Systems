function bc=ac6bcfx_jac(p,u) 
% BC for geometry with six edges, mixed Neumann and Dirichlet, x-dependent 
sf=10^4; % stiff spring factor for Dirichlet part
qdj = mat2str(sf);gdj ='0'; % pseudo-Dirichlet: 
qn = '0'; gn = '0'; % Neumann part
% ordering in geometry from rec: bottom, right, top, left
bc = gnbcs(p.nc.neq,qn,gn,qn,gn,qn,gn,qn,gn,qdj,gdj,qn,gn);