function [x,p]=lssgmres(A,b,p) % gmres interface for SH1D, see afun 
n=size(b,1); try; ittol=p.ittol; limax=p.limax; catch; ittol=1e-8; limax=n; end 
maxit=min(limax,n); L=p.mat.prec; % preconditioner (diagonal, shifted multipl)
tic;[x,flag,relres,iter]=gmres(@afun,b,[],ittol,maxit,L,L');iter=iter(2);t1=toc; 
if p.sw.verb>2; fprintf('gmres-flag=%i, relres=%g, i=%i, time=%g\n',...
                        flag,relres,iter,t1); end  
% following line is a typical 'fallback' (but not here, cause A is not set)
if flag>0; fprintf('gmres failed, using lss\n'); tic; x=A\b; toc, end 