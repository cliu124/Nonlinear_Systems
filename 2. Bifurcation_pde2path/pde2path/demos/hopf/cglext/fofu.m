function f=fofu(p,u) % forcing function, using: 
% p.T=period, p.t=current time (normalized to (0,1)), set in hosrhs.m 
% with error catching for steady problem (or if forcing is not set) 
par=u(p.nu+1:end); 
try; pa=par(8); T=p.T; t=T*p.t; tc=par(9); catch; pa=0; T=0; t=0; tc=0; end 
x=getpte(p); n=p.nu/2; x=x(1:n)'; f=pa*tanh(10*(t-tc*T))*sin(x); 