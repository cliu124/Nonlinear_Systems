function [y,T,lam,res,iter,A,p,jac]=honloop_np(p,y1,T1,lam1,ds)
% honloop_np: Newton loop for Hopf in natural parametrization (fixed lam) 
% assemble via TOM, but solve here, either with fsolve or Newton-loop 
if p.fsol.fsol==1; 
    [y,T,lam,res,iter,A,p,jac]=fsol_np(p,y1,T1,lam1,ds); 
    return; 
end
iter=0; y=y1; T=T1; lam=lam1; p=setlam(p,lam); % set lam to predictor
nu=p.nu; na=nu; tl=p.hopf.tl; tol=p.nc.tol; 
par=getaux(p); s=par(p.hopf.ilam); % active aux vars
alpha=1; almin=alpha/4; % minimal damping 
[f,jac,f_T,f_lam,f_s]=tomassempbc(p,y,T,lam); 
[pc,pc_y]=hopc(p,y,T,lam); %arcl=hoarcl(p,y,T,lam,ds);
qf=p.hopf.qfh(p,y); qhder=getqhder(p,y); 
F=[f; pc; qf]; 
fprintf('initial (||f||, pc,qf,a)=(%g, %g, %g ,%g)\n',norm(f,p.sw.norm),pc,qf,s);
res0=norm(F,p.sw.norm); res=res0; % pause  % the residual 
ytom=reshape(y1(1:na*tl),na,tl); % convert to tom-input 
A=gethoA_np(p,jac,f_T,pc_y,f_s,qhder); 
if(res<tol); return; end; % initial res. small, do nothing! 
% now start the actual loop
stepok=1; % stepok=1 indicates that step of size alpha decreased residual  
while(abs(res0)>tol && iter<p.nc.imax && stepok) 
  [upd,p]=p.fuha.blss(A,F,p); 
  stepok=0; yv=reshape(y(1:na,:),na*tl,1); % make y-vector
  while(stepok==0 && alpha>=almin)   
  y1=yv-alpha*upd(1:na*tl); 
  T1=T-alpha*upd(na*tl+1); 
  s1=s-alpha*upd(na*tl+2);% a,a1
  p.u(p.nu+p.hopf.ilam)=s1; 
  ytom=reshape(y1(1:na*tl),na,tl); % convert to tom-input 
  f=tomassempbc(p,ytom,T1,lam1);  % setlam(p,lam1) done in tomassempbc! 
  [pc,pc_y]=hopc(p,ytom,T1,lam1); % arcl=hoarcl(p,ytom,T1,lam1,ds); 
  qf=p.hopf.qfh(p,y); qfj=p.hopf.qfhder(p,y); F=[f; pc; qf];
  fprintf('al=%g, (||f||, pc, qf,a)=(%g,%g,%g,%g)\n',alpha,norm(f,p.sw.norm),pc,qf,s1);
  res=norm(F,p.sw.norm); if p.sw.verb>1; disp(['res=' mat2str(res,6)]); end % check new tom-resi 
  %pause 
  if(res<res0); % good step  
    if(res<res0/2 && alpha<1) alpha=alpha*2; end % very good step  
    stepok=1; y=ytom; T=T1; lam=lam1; s=s1; res0=res;
    [f,jac,f_T,f_lam, f_s]=tomassempbc(p,y,T,lam); 
    [pc,pc_y]=hopc(p,y,T,lam); qhder=getqhder(p,y); 
    A=gethoA_np(p,jac,f_T,pc_y,f_s,qhder); 
  else alpha=alpha/2; % bad step, try smaller alpha
  end % res<res0
  end % while stepok==0
  iter=iter+1;
end % while res<p.nc.tol && stepok 
%fprintf('al=%g, (||f||, pc, arcl)=(%g,%g,%g)\n',alpha,norm(f,p.sw.norm),pc,arcl);
% ----------------------------------------------------------postprocessing
if(p.sw.verb>1); if(alpha<1) fprintf('\nhonloopext: damp alpha=%g, res=%g, ds=%g\n',...
            alpha,res,ds); end; end; % inform user if damping was used 
res=res0; y=ytom; lam=lam1; 
