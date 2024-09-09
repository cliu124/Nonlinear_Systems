function p=oosetfemops(p) % for Brusselator 
grid=p.pdeo.grid; par=p.u(p.nu+1:end); Du=par(5); Dv=par(6); Dw=par(7); 
[K,M,~]=p.pdeo.fem.assema(grid,1,1,1); 
N=sparse(grid.nPoints, grid.nPoints); 
p.mat.K=[[Du*K N N];[N Dv*K N]; [N N Dw*K]]; 
p.mat.M=kron([[1,0,0];[0,1,0]; [0,0,1]],M); 
end 