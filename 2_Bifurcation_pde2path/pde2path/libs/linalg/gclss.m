function [x,p]=gclss(A,b,p) 
% gclss: customized LinearSystemSolver for global coupling. 
% To solve: Av-c(u)<h_u(u)*v>=b 
% With <u>=av(u), write (A-cvec*avec)x=b and solve by Sherman-Morrison
global p2pglob; % cvec and avec set in sG and sGjac, resp 
cvec=p2pglob.cvec; avec=p2pglob.avec; 
y=A\b; z=A\cvec; al0=avec*z; al=avec*y/(1-al0); x=y+al*z; 
return % comment out to check Sherman-Morrison formula
c1=(A-Mc*avec)*x; c2=(A-Mc*avec)*y; 
err1=norm(c1-b); err2=norm(c2-b); 
fprintf('lss: |b|=%g, |x|=%g, |y|=%g, SM-err=%g, noSM-err=%g\n',...
    norm(b),norm(x),norm(y),err1,err2); pause 