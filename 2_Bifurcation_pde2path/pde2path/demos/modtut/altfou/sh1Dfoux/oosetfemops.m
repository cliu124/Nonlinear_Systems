function p=oosetfemops(p) % SH via dct, u in x-space, diff matrix L=F'*mu*F  
kfx=pi/p.lx; n=p.np; p.mat.M=speye(n); 
kvx=kfx*[0:n-1]'; F=dctmtx(n); % wave-nr and dct-matrix 
dd=spdiags((1-kvx.^2).^2,0,n,n); % multipliers of L  
p.mat.L=F'*dd*F; 
