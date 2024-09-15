function p=oosetfemops(p) % for pollution demo
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
p.mat.K=kron(diag([p.d1 p.d2 -p.d1 -p.d2]),K); % note the anti-diff in 3 and 4
p.mat.M=kron(diag([1 1 1 1]),M); 