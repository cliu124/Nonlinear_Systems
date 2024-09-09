function p=oosetfemops(p)  % ac3D 
gr=p.pdeo.grid; p.nc.sf=1e3; % stiff spring constant for DBC via Robin-BC
% see sGws for param dependent BCs! 
[p.mat.K,p.mat.M,~]=p.pdeo.fem.assema(gr,1,1,1); % indep. of BC 