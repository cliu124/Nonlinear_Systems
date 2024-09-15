function [qL,qU,qD]=stanqM(p)
% stanqM: standard extension (by zero borders) of M for q (constraints) 
qL=zeros(p.nc.nq,p.nu); qU=zeros(p.nu,p.nc.nq); qD=0*speye(p.nc.nq);  
