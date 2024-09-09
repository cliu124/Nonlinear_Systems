function r=nleq1resi(p,au)
% NLEQ1RESI: resi wrapper for NLEQ1
u=au2u(p,au); r=resi(p,u); 
