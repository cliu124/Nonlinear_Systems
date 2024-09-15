function q=qfh(p,y) % aux eqns in Hopf, here: sum up shifts wrt u0
n=p.nu; tl=size(p.hopf.y,2); q=0;  % phase constr, useful to define 'on average' 
u0x=p.u0x(1:p.nu); 
for i=1:tl; u=y(1:n,i); q=q+u0x'*u; end 