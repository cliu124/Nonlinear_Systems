function y=myeigsfu(p,A,B,sig,b) 
% gcafun: solver to be called in eigs; calls ilueigs 
A=A-sig*B; p.pcgtol=1e-8; y=lssgmres(A,b,p); p.pcgtol=1e-8; 