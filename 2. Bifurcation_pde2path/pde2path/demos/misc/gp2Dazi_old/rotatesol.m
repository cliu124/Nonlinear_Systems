function p=rotatesol(p,al) % rotate sol and check residual
r=resi(p,p.u); res1=norm(r);
u1=p.u(1:p.np); u2=p.u(p.np+1:2*p.np); par=p.u(p.nu+1:end); 
v1=cos(al)*u1-sin(al)*u2; v2=sin(al)*u1+cos(al)*u2;
p.u=[v1;v2;par];r=resi(p,p.u); res2=norm(r);
fprintf('Phase rotation of sol by alpha=%g, old resi=%g, new resi=%g\n',al,res1,res2); 