function p=oosetfemops(p) % SH by dct, matrix free, multipl. stored in p2pglob
n=p.np; p.mat.M=speye(n); kfx=pi/p.lx; kvx=kfx*[0:n-1]'; % M and wave numbers 
global p2pglob; p2pglob.mu=(1-kvx.^2).^2; % multiplier, used in sG and afun 
p.mat.prec=spdiags(sqrt((1-kvx.^2).^2+1),0,n,n); % prec, used in lssgmres 