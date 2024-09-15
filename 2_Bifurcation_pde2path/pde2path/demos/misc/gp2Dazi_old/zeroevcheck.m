function zeroevcheck(p) % calc. "zero" eigenvector of G_u for GP
% and check that it's [-v u]
r=[]; Gu=getGu(p,p.u,r); cest=condest(Gu); 
vs=size(Gu,1); p.evopts.v0=ones(vs,1)/vs;
[ev,mu]=eigs(Gu,1,0,p.evopts); 
u1=p.u(1:p.np); u2=p.u(p.np+1:2*p.np); u1=u1/norm(u1); u2=u2/norm(u2);
v1=ev(1:p.np); v2=ev(p.np+1:2*p.np); v1=v1/norm(v1); v2=v2/norm(v2);
n1=norm(u1+v2); n2=norm(u2-v1); n3=norm(u1-v2); n4=norm(u2+v1); 
% print the minimum since we don't know which direction of ev is chosen
fprintf('zero-EVal mu=%g, cond-est=%g, EVec-check: %g %g\n',...
    mu,cest,min(n1,n3),min(n2,n4));