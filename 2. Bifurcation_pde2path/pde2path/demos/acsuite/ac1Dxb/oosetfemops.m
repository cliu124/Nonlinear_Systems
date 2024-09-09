function p=oosetfemops(p) % for AC1Dxb: x-dependent diffusion in sG and sGjac,
gr=p.pdeo.grid; [K,M,~]=p.pdeo.fem.assema(gr,1,1,1); p.mat.M=M; p.mat.K=K; 
p.mat.Kx=p.pdeo.fem.convection(gr,1); % convection matrix 