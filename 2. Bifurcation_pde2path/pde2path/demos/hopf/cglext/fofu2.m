function f=fofu2(p,u) % forcing function, using: 
% p.T=period (as set in hosrhs.m) and 
% p.t=current time (normalized to (0,1))
% with error catching for steady problem (or if forcing pars are not set) 
par=u(p.nu+1:end); c5=par(5); 
try; t=p.T*p.t; catch; t=0; end 
x=getpte(p); x=x'; uv=cos(par(7)*x);
f=c5*cos(t)*uv;