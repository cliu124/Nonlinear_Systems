function [f,jac,f_T,f_lam,f_a]=huassempbc(p,y,T,lam)
% huassempbc: assemble rhs and jac for hopf with pBC, pad y accordingly
%
% [f,jac,f_T,f_lam]=huassempbc(p,ya,T,lam) 
switch p.hopf.dsw; 
    case 1; afu=@huassem1; % forward FD
    case 2; afu=@huassem2; % central FD 
    case 3; afu=@huassem3; % backward FD
    case 4; afu=@huassem4; % TOM
    case 5; afu=@huassem5; % impl mid-point
    case 6; afu=@huassem6; % impl mid-point for G, on point for q
    case 10; afu=@huassem10; % user 
end
try; freeT=p.hopf.freeT; catch; freeT=1; end % check if T is free (default)
nu=p.nu; nq=p.nc.nq; na=nu+nq; tl=p.hopf.tl; if ~isfield(p.hopf,'jac'); p.hopf.jac=1; end 
p=setlam(p,lam); par=p.u(nu+1:end); 
if nargout==1; f=afu(p,y,T,par);  return; end  %  **** if only f is required ****
[f,jac]=afu(p,y,T,par); mclf(4); plot(f(1:3*nu),'-*'); title('f'); 
del=p.nc.del; % FD approximate f_T 
f1=afu(p,y,T+del,par); f_T=(f1-f)/del; %mclf(5); 11, plot(f_T(1:3*nu+10),'-*'); title('f_T'); 
p=setlam(p,lam+del); par=p.u(nu+1:end); % FD approximate \pa_lam f 
f1=afu(p,y,T,par);   
f_lam=(f1-f)/del; mclf(6); % plot(f_lam(1:nu)); pause 
p=setlam(p,lam); 
f_a=0; % Hopf aux vars
try; nqh=length(p.hopf.ilam); catch; nqh=0; end; 
if nqh>0  % aux Hopf-eqns-derivative 
  f_a=zeros(tl*na,nqh); 
  for i=1:nqh  
    ps=par(p.hopf.ilam(i)); par(p.hopf.ilam(i))=ps+del;
    f1=afu(p,y,T,par);   
    f_a(:,i)=(f1-f)/del; 
    par(p.hopf.ilam(i))=ps;  
  end
end 
