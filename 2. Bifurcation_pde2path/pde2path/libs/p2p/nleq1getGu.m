function Gu=nleq1getGu(p,au,r)
% NLEQ1GETGU: wrapper for getGU 
u=au2u(p,au); Gu=getGu(p,u,r); 
