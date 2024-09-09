function [x,p]=lssgmres(A,b,p) % gmres for sh1D  
n=size(A,1); maxit=min(100,n); 
L=p.mat.prec; try; ittol=p.ittol; catch;  ittol=1e-8; end 
tic; [x,flag,relres,iter]=gmres(A,b,[],ittol,maxit,L,L'); iter=iter(2); t1=toc; 
% tic; x=A\b; toc % just comp., gmres wins for SH for nu>5000
if p.sw.verb>2; fprintf('gmres-flag=%i, relres=%g, i=%i, time=%g\n',flag,relres,iter,t1); end  
if flag>0; fprintf('gmres failed, using lss\n'); tic; x=A\b; toc, end 
