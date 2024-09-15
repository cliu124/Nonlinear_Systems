function A=gethoA(p,jac,f_T,f_lam,pc_y,tau,f_a, qhder)
% gethoA: build extended Jac for Hopf-arclength setting
nqh=p.hopf.nqh; try; nqh2=length(p.hopf.ilam); catch; nqh2=0; end 
try; freeT=p.hopf.freeT; catch; freeT=1; end % check if T is free (default) 
try; pcsw=p.hopf.pc; catch; pcsw=1; end % pcsw=1: with phase cond. 
xi=p.hopf.xi; tw=p.hopf.tw; nu=p.nu; na=nu+p.nc.nq; tl=p.hopf.tl; 
if freeT; % and for now pc=1 
 pc_T=0; 
 A=[[jac, f_T, f_lam];
    [pc_y, pc_T, 0];  
    [xi*tau(1:na*tl), (1-xi)*tw*tau(na*tl+1), (1-xi)*(1-tw)*tau(na*tl+2)]];
 if nqh>0;  A=[[A [f_a; zeros(2,nqh2)]]; qhder]; end 
else % T fixed, no \pa_T(G) column in A 
 if pcsw; % with phase-cond. 
  A=[[jac, f_lam, f_a];
    [pc_y, 0, zeros(1,nqh2)];  
    [xi*tau(1:na*tl), (1-xi)*tau(na*tl+1:na*tl+1+nqh2)]];
    if nqh>0; A=[A; qhder]; end    
 else % no pc 
  A=[[jac, f_lam];
    [xi*tau(1:na*tl), (1-xi)*tau(na*tl+1)]];  
    if nqh>0; A=[A [f_a; zeros(1,nqh2)]; qhder];  end      
 end
end