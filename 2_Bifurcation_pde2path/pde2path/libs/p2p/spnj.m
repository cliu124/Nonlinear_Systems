function r=spnj(t,u)
% spnj: interface routine to compute pa_u(G_u\phi) via numjac 
% r=spnj(t,u)
global pj phij; 
p=pj; u1=[u; p.u(2*p.nu+1:end)]; % append pars
r1=pderesi(pj,u1); Gu=getGupde(p,u1,r1); r=Gu*phij; 