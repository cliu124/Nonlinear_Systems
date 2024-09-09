function [x,p]=lssgmres(A,b,p) % matrix free LSS via afun, A is dummy here 
n=size(A,1); try; ittol=p.ittol; limax=p.limax; catch; ittol=1e-8; limax=n; end 
maxit=min(p.limax,n); L=p.mat.prec; 
tic; [x,flag,relres,iter]=gmres(@afun,b,[],ittol,maxit,L,L'); iter=iter(2); t1=toc; 
if p.sw.verb>2; fprintf('gmres-flag=%i, relres=%g, i=%i, time=%g\n',...
                        flag,relres,iter,t1); end  
if flag>0; fprintf('gmres failed, using lss\n'); tic; x=A\b; toc, end 
