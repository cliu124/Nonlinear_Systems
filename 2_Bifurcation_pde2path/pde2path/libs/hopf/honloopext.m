function [y,T,lam,a,res,iter,A,p,jac]=honloopext(p,y1,T1,lam1,a1,ds,varargin)
% honloopext: Newton loop for Hopf in arclength
%
% updated for cases 
% without t-transl. phase-condition: p.hopf.pc=0  (default 1) 
try; freeT=p.hopf.freeT; catch; freeT=1; end % check if T is free (default) 
try; pcsw=p.hopf.pc; catch; pcsw=1; end  % check if t-PC is used (default) 
try; nqh=p.hopf.nqh; catch; nqh=0; end 
try; nqh2=length(p.hopf.ilam); catch; nqh2=0; end 
try; a=p.u(p.nu+p.hopf.ilam); catch; a=0; end; 
np=p.np; 
iter=0; y=y1; T=T1; lam=lam1; tol=p.nc.tol; 
p=setlam(p,lam); % set lam to predictor
nu=p.nu; tl=p.hopf.tl; 
alpha=1; try; almin=p.nc.almin; catch;  almin=0.5; end % minimal damping 
%[f,jac,f_T,f_lam,f_a]=tomassempbc(p,y,T,lam); 
[f,jac,f_T,f_lam,f_a]=hoassempbc(p,y,T,lam); 
pc=[]; pc_y=[]; if pcsw; [pc,pc_y]=hopc(p,y,T,lam); end 
qf=[]; qhder=[]; if nqh>0; qf=p.hopf.qfh(p,y); qhder=getqhder(p,y); end; %qf, pause 
arcl=hoarcl(p,y,T,lam,ds,a);  %size(f), size(pc), size(arcl), size(qf), pause 
F=[f; pc; arcl; qf]; res0=norm(F,p.sw.norm); res=res0; 
if p.sw.verb>1;  
 if pcsw; fprintf('initial (||f||, pc, arcl)=(%g,%g,%g)\n',norm(f,p.sw.norm),pc,arcl); 
 else; fprintf('initial (||f||, arcl)=(%g,%g)\n',norm(f,p.sw.norm),arcl); end
 if nqh>0; fprintf('qh=%g\n',qf); end 
end
%size(jac), size(qhder), size(p.hopf.tau), pause 
A=gethoA(p,jac,f_T,f_lam,pc_y,p.hopf.tau,f_a,qhder); %mclf(8); spy(A), %pause 
if(res<tol); return; end; % initial res. small, do nothing! 
% now start the actual loop
stepok=1; % stepok=1 indicates that step of size alpha decreased residual  
while(abs(res0)>tol && iter<p.nc.imax && stepok); % size(F), size(A)
  [upd,p]=p.fuha.blss(A,F,p); %upd(nu*tl+1) 
  stepok=0; yv=reshape(y(1:nu,:),nu*tl,1); % make y-vector
  while(stepok==0 && alpha>=almin)   
  y1=yv-alpha*upd(1:nu*tl); 
  if freeT; tupd=upd(nu*tl+1); T1=T-alpha*tupd; % tupd, pause 
     lam1=lam-alpha*upd(nu*tl+2); 
     if nqh2>0; a=p.u(p.nu+p.hopf.ilam);  
         a1=a-alpha*upd(nu*tl+3:end); 
         p.u(p.nu+p.hopf.ilam)=a1; 
     end 
  else; T1=T;  lam1=lam-alpha*upd(nu*tl+1); %nqh2, pause 
    if nqh2>0;  a=p.u(p.nu+p.hopf.ilam); 
        a1=a-alpha*upd(nu*tl+2:end);   p.u(p.nu+p.hopf.ilam)=a1; 
    end
  end 
  ytom=reshape(y1(1:nu*tl),nu,tl); % convert to tom-input 
  f=hoassempbc(p,ytom,T1,lam1);  % setlam(p,lam1) done in tomassempbc! 
 % f(3*np+1), pause 
  if pcsw; [pc,pc_y]=hopc(p,ytom,T1,lam1); end 
  if freeT; 
     arcl=hoarcl(p,ytom,T1,lam1,ds,a1); 
     if nqh>0; qf=p.hopf.qfh(p,ytom); end
  else 
   if nqh2>0; p.u(p.nu+p.hopf.ilam)=a; arcl=hoarcl(p,ytom,T1,lam1,ds,a1); 
       p.u(p.nu+p.hopf.ilam)=a1; 
   else; arcl=hoarcl(p,ytom,T1,lam1,ds,a1); 
   end
  end 
  F=[f; pc; arcl]; if nqh>0; F=[F; qf]; end 
  if p.sw.verb>1; 
    %  alpha, size(qf), pause 
    if nqh==0; 
     if pcsw; fprintf('al=%g, (||f||, pc, arcl)=(%g,%g,%g)\n',alpha,norm(f,p.sw.norm),pc,arcl);
     else; fprintf('al=%g,  (||f||, arcl)=(%g,%g)\n',alpha,norm(f,p.sw.norm),arcl); 
     end
    else
     if pcsw; fprintf('al=%g, (||f||, pc, arcl,qh)=(%g,%g,%g,%g)\n',alpha,norm(f,p.sw.norm),...
             pc,arcl,norm(qf,'inf'));
     else; fprintf('al=%g,  (||f||, arcl, qh)=(%g,%g,%g)\n',alpha,norm(f,p.sw.norm),arcl,qf); 
     end 
    end
    %if nqh>1; fprintf('qh=%g\n',qf);  end 
  end
  res=norm(F,p.sw.norm); % check new tom-resi 
  if(res<res0); % good step       
    if(res<res0/2 && alpha<1) alpha=alpha*2; end % very good step  
    stepok=1; y=ytom; T=T1; lam=lam1; res0=res;
    [f,jac,f_T,f_lam,f_a]=hoassempbc(p,y,T,lam); 
    if pcsw; [pc,pc_y]=hopc(p,y,T,lam); end 
    if nqh==0; A=gethoA(p,jac,f_T,f_lam,pc_y,p.hopf.tau,f_a); 
    else; qhder=getqhder(p,y); 
     A=gethoA(p,jac,f_T,f_lam,pc_y,p.hopf.tau,f_a,qhder); 
    end
  else alpha=alpha/2; % bad step, try smaller alpha
  end % res<res0
  end % while stepok==0
  iter=iter+1;
end % while res<p.nc.tol && stepok 
% ----------------------------------------------------------postprocessing
if(p.sw.verb>1); if(alpha<1) fprintf('\nhonloopext: damp alpha=%g, res=%g, ds=%g\n',...
            alpha,res,ds); end; end; % inform user if damping was used 
res=res0; try; y=ytom; catch; y=y1; end 

