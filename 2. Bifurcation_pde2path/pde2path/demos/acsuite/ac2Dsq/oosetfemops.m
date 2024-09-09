function p=oosetfemops(p) % ac2D 
gr=p.pdeo.grid; 
[p.mat.K,p.mat.M,~]=p.pdeo.fem.assema(gr,1,1,1);  
