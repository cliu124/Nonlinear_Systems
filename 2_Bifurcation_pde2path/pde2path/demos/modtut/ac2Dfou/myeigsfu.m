function y=myeigsfu(p,A,B,sig,b) 
% function to be called in eigs; 
A=A-sig*B; y=lssgmres(A,b,p); 