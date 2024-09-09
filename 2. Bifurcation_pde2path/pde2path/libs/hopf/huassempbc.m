function [f,jac,f_T,f_lam,f_a]=huassempbc(p,y,T,lam)
% huassempbc: assemble rhs and jac for hopf with pBC, pad y accordingly
%
% [f,jac,f_T,f_lam]=huassempbc(p,ya,T,lam) 
try; freeT=p.hopf.freeT; catch; freeT=1; end % check if T is free (default)
nu=p.nu; nq=p.nc.nq; na=nu+nq; tl=p.hopf.tl; if ~isfield(p.hopf,'jac'); p.hopf.jac=1; end 
p=setlam(p,lam); par=p.u(nu+1:end); t=p.hopf.t; 
if nargout==1; f=huassem(p,T,y,par); return; end  %  **** if only f is required ****
[f,jac]=huassem(p,T,y,par); % mclf(4); plot(f(1:2*nu)); pause
del=p.nc.del; % FD approximate f_T 
f1=huassem(p,T+del,y,par); f_T=(f1-f)/del; 
p=setlam(p,lam+del); par=p.u(nu+1:end); % FD approximate \pa_lam f 
f1=huassem(p,T,y,par);   
f_lam=(f1-f)/del; %mclf(4); plot(f_lam(1:4*nu)); pause 
p=setlam(p,lam); 
f_a=0; % Hopf aux vars
try; nqh=length(p.hopf.ilam); catch; nqh=0; end; 
if nqh>0  % aux Hopf-eqns-derivative 
  f_a=zeros(tl*na,nqh); 
  for i=1:nqh  
    ps=par(p.hopf.ilam(i)); par(p.hopf.ilam(i))=ps+del;
    f1=huassem(p,T,y,par);   
    f_a(:,i)=(f1-f)/del; 
    par(p.hopf.ilam(i))=ps;  
  end
end 
