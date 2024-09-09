function p=oosetfemops(p) % SH via dct matrix F and multipl. matrix L 
n=p.np; p.mat.M=speye(n); % mass matrix is Identity 
kfx=pi/p.lx; kvx=kfx*[0:n-1]'; p.mat.F=dctmtx(n); % wave-nr and dct-matrix 
p.mat.L=spdiags((1-kvx.^2).^2,0,n,n); % multipliers as a diag matrix 

