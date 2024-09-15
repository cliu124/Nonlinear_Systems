function y=gcafun(p,A,B,sig,b) 
% gcafun: solver to be called in eigs; calls SM solver gclsseigs
A=A-sig*B; y=gclsseigs(A,b,p); 