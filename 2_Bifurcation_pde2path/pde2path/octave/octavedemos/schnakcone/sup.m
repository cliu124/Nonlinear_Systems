close all; clear classes
%% Schnakenberg on pointed cone (a=0.2), scaled by 1/eps 
dir='h1'; nx=50; ell=1.2; eps=0.04; lam=3.225; sig=0; D=60; a=0.2; s=0; dsmin=0.005; nx=20; 
p=[]; par=[lam sig D a s eps]; p=schnakcinit(p,nx,par,ell); p.sol.ds=-dsmin;  
p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.lammin=3; p.nc.usrlam=[3 3.1 3.2]; p0=p; 