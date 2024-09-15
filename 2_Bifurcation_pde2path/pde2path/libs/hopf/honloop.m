function [y,T,lam,res,iter,A,p]=honloop(p,y1,T1,lam1)
% honloop: Newton loop for Hopf at fixed lam, for now nqh=0
try; pcsw=p.hopf.pc; catch; pcsw=1; end 
try; freeT=p.hopf.freeT; catch; freeT=1; end % check if T is free (default) 
%if ~pcsw; [y,T,lam,res,iter,A,p]=nahonloop(p,y1,T1,lam1); return; end  
iter=0; y=y1; T=T1; lam=lam1; p=setlam(p,lam); nu=p.nu; tl=p.hopf.tl; 
%mbws=p.nc.mbw; p.nc.mbw=1; % save and reset border width
alpha=1; try; almin=p.nc.almin; catch;  almin=0.5; end % minimal damping 
[f,jac,f_T,f_lam]=tomassempbc(p,y,T,lam); 
pc=[]; if pcsw; [pc,pc_y]=hopc(p,y,T,lam); end 
F=[f; pc];  
if p.sw.verb>1;  
 if pcsw; fprintf('initial (||f||, pc)=(%g,%g)\n',norm(f,p.sw.norm),pc); 
 else; fprintf('initial ||f||=%g\n',norm(f,p.sw.norm)); end
end
res0=norm(F,p.sw.norm); res=res0;  % the residual 
if pcsw; A=[[jac, f_T]; [pc_y, 0]] ; else A=jac; end 
tol=p.nc.tol; 
if(res<tol); return; end; % initial res. small, do nothing! 
% now start the actual loop
stepok=1; % stepok=1 indicates that step of size alpha decreased residual! 
while(abs(res0)>tol && iter<p.nc.imax && stepok)     
  [upd,p]=p.fuha.lss(A,F,p); % (for lssbel has bw=bw of big_bel-1) 
  stepok=0; yv=reshape(y(1:nu,:),nu*tl,1); % make y-vector
  while(stepok==0 && alpha>=almin)   
  y1=yv-alpha*upd(1:nu*tl);
  if freeT; T1=T-alpha*upd(nu*tl+1); else; T1=T; end 
  ytom=reshape(y1(1:nu*tl),nu,tl); % convert to tom-input 
  f=tomassempbc(p,ytom,T1,lam); 
  if pcsw; [pc,pc_y]=hopc(p,ytom,T1,lam); end 
  F=[f; pc];  res=norm(F,p.sw.norm); % check new tom-resi
  if p.sw.verb>1;  
    if pcsw; fprintf('al=%g, (||f||, pc)=(%g,%g)\n',alpha,norm(f,p.sw.norm),pc); 
    else; fprintf('al=%g,  ||f||=%g\n',alpha,norm(f,p.sw.norm)); end
  end  
  if(res<res0); % good step       
    if(res<res0/2 && alpha<1) alpha=alpha*2; end % very good step  
    stepok=1; y=ytom; T=T1; res0=res;
    [f,jac,f_T,f_lam]=tomassempbc(p,y,T,lam); 
    if pcsw; [pc,pc_y]=hopc(p,y,T,lam); end 
    if pcsw; A=[[jac, f_T]; [pc_y, 0]]; else; A=jac; end 
  else alpha=alpha/2; % bad step, try smaller alpha
  end % res<res0
  end % while stepok==0
  iter=iter+1;
end % while res<p.nc.tol && stepok 
%fprintf('al=%g, (||f||, pc, arcl)=(%g,%g,%g)\n',alpha,norm(f,p.sw.norm),pc,arcl);
% ----------------------------------------------------------postprocessing
if(p.sw.verb>1); if(alpha<1) fprintf('\nhonloop: damp alpha=%g, res=%g,\n',...
            alpha,res0); end; end; % inform user if damping was used 
res=res0; try; y=ytom; catch; end; %p.nc.mbw=mbws; 
%fprintf('honloop, lam=%g, res=%g\n', lam,res); 
