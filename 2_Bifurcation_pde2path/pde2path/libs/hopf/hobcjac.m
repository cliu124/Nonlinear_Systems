function [ja,jb]=hobcjac(ya,yb,p,par) 
% hobc: BC jacobian for Hopf, for tomsol via tom 
nu=p.nu; sd=spdiags(ones(nu,1),0,nu,nu); 
ja=[[sd -sd zeros(nu,1)]; 
    [0*sd 0*sd zeros(nu,1)]
    [p.hopf.u0dot' zeros(1,nu+1)]]; 
jb=[[0*sd 0*sd zeros(nu,1)]; 
    [sd -sd zeros(nu,1)]; 
     zeros(1,2*nu+1)]; 



