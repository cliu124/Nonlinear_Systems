function y=gcafunilu(p,A,B,sig,b) 
% gcafun: solver to be called in eigs; calls ilueigs 
A=A-sig*B; y=ilueigs(A,b,p); 