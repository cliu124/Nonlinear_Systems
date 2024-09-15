function p=oosetfemops(p) % FFT based, but u in x-space 
kf=pi/p.lx; n=p.np; p.mat.M=speye(n); 
kv=kf*[0:n-1]'; % F-vectors (normalized to (0,lx)) 
p.mat.k2=kv.^2;  % F multipliers 
F=dctmtx(n);    % discrete-cos-transf. matrix 
p.mat.L=F'*spdiags(p.mat.k2,0,n,n)*F;