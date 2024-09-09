function [x,p]=lssgmres(A,b,p)
n=size(A,1); maxit=min(100,n); 
L=p.mat.prec; try; ittol=p.ittol; catch;  ittol=1e-8; end 
tic; 
[x,flag,relres,iter]=gmres(A,b,[],ittol,maxit,L,L');
t1=toc; 
if p.sw.verb>2; fprintf('gmres-flag=%i, relres=%g, i1=%i, i2=%i, time=%g\n',flag,relres,iter(1),iter(2),t1); end  
if flag>0; fprintf('gmres failed, using lss\n'); tic; x=A\b; toc, end 
