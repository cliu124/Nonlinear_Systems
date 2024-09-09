function [ja,jb]=hodummybcjac(ya,yb,p,par,T)  
% hodummybcjac: dummy BC jac for Hopf, used (but irrelevant) in tomassempbc
nd=p.nu+p.nc.nq; sd=spdiags(ones(nd,1),0,nd,nd); 
ja=sd; jb=sd; 


