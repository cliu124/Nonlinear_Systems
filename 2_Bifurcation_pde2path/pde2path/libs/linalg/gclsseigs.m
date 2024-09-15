function [x,p]=gclsseigs(A,b,p) 
% gclsseigs: customized LinearSystemSolver for GC, eigs-version  
global p2pglob; cvec=p2pglob.cvec; avec=p2pglob.avec;
[y,p]=lsslueigs(A,b,p); [z,p]=lsslueigs(A,cvec,p); 
%[y,p]=lss(A,b,p); [z,p]=lss(A,cvec,p); 
al0=avec*z; al=avec*y/(1-al0); x=y+al*z; 
return % comment out to check Sherman-Morrison formula
c1=(A-lam*cvec*avec)*x; c2=(A-lam*cvec*avec)*y; 
err1=norm(c1-b); err2=norm(c2-b); 
fprintf('lss: |b|=%g, |x|=%g, |y|=%g, SM-err=%g, noSM-err=%g\n',...
    norm(b),norm(x),norm(y),err1,err2); pause 