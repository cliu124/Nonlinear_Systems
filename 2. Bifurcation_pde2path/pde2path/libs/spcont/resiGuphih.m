function r=resiGuphih(t,au)
% resiGuphi: interface for numjac for \pa_u(Gu*phi) 
% currently only works with p.sw.jac=1; 
global pGu % the parameters and psi from pGu
p=pGu;
ua=[au; p.u(p.nu+1:end)]; phir=p.phir; phii=p.phii;
r0=pderesi(p,ua); Gu=getGupde(p,ua,r0); 
r=[Gu*phir;Gu*phii];  
