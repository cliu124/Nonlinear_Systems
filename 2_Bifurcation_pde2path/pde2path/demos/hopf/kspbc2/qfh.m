function q=qfh(p,y) % aux eqns in Hopf, here: sum up shifts wrt u0
par=p.u(p.nu+1:end); m=par(2); n=p.nu/2; 
q1=sum(p.mat.Ms*y(1:n,1))/p.vol-m; % mass constraint (at initial slice) 
tl=size(p.hopf.y,2); q2=0;  % phase constr, useful to define 'on average' 
if isfield(p,'u0x'); u0x=p.u0x(1:n); 
else u0=p.u(1:n); u0x=p.mat.Ms*(p.mat.Kx*u0); end
for i=1:tl; u=y(1:n,i); q2=q2+u0x'*u(1:n); end 
q=[q1;q2]; 