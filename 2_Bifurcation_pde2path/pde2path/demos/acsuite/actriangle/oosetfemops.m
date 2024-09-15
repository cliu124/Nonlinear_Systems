function p=oosetfemops(p) % ac2D 
gr=p.pdeo.grid; p.nc.sf=1e4; % stiff spring constant for DBC via Robin-BC
[p.mat.K,p.mat.M,~]=p.pdeo.fem.assema(gr,1,1,1);  % indep. of BC 
