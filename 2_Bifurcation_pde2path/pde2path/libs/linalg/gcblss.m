function [x,p]=gcblss(A,b,p) 
% : customized blss for GC  
global p2pglob; cvec=p2pglob.cvec; avec=p2pglob.avec;
ecvec=[cvec; 0]; eavvec=[avec 0]; % extented global vectors 
y=A\b; z=A\ecvec; al0=eavvec*z; al=eavvec*y/(1-al0); x=y+al*z; 
return % comment out to check Sherman 
c1=(A-ecvec*eavvec)*x; c2=(A-ecvec*eavvec)*y; 
err1=norm(c1-b); err2=norm(c2-b); 
fprintf('blss: SM-err=%g, noSM-err=%g\n', err1,err2);  pause