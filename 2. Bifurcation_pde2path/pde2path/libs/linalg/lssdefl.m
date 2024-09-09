function [x,p]=lssdefl(A,b,p) 
% lssdefl: customized LinearSystemSolver for deflation 
% To solve: Av-(c*a)*v=b 
global p2pglob; % cvec and avec set in sG and sGjac, resp 
nu=p.nu; cvec=p2pglob.cvec(1:nu); avec=p2pglob.avec(1:nu); 
y=A\b; 
if norm(avec,2)>1e-8; z=A\cvec; al0=avec*z; al=avec*y/(1-al0); x=y+al*z; %al
else x=y; 
end
return % comment out to check Sherman-Morrison formula
c1=(A-Mc*avec)*x; c2=(A-Mc*avec)*y; 
err1=norm(c1-b); err2=norm(c2-b); 
fprintf('lss: |b|=%g, |x|=%g, |y|=%g, SM-err=%g, noSM-err=%g\n',...
    norm(b),norm(x),norm(y),err1,err2); pause 