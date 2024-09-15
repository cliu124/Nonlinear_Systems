function p=oosetfemops(p) % FFT based, u in Fourier-space 
kf=pi/p.lx; n=p.np; p.mat.M=speye(n); 
kv=kf*[0:n-1]'; % F-vectors (normalized to (0,lx)) 
p.mat.k2=kv.^2; p.mat.F=dctmtx(n); % store multipliers,  and dct matrix 
