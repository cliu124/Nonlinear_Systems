function r=resiGuTpsi(t,au)
% resiGuTpsi: interface for numjac for \pa_u(Gu'*psi) 
global pGu % the parameters and psi from pGu
p=pGu;
ua=[au; p.u(p.nu+1:end)]; psipde=p.psipde;
r0=pderesi(p,ua); Gu=getGupde(p,ua,r0); r=Gu'*psipde; 
%13, norm(r-pGu.r1), pause 
%11, ua(p.nu+1:end)', ua(1:4)', r0(1:4)', r(1:4)', pause 
