function [u,res,iter,Gu,Glam,p]=nloop(p,u1)
% NLOOP: Newton loop with damping according to settings, 
%  see in particular stanparam. 
%
%  [u,res,iter,Gu,Glam,p]=nloop(p,u1)
%
% * u1 is initial guess 
% * Returns resulting u, res=residuum, iter=number of iterations, Gu,Glam
%   derivatives of G. 
% * Note: returns full vector u (including auxiliary variables)! 
% * Important settings: p.fsol.fsol, p.nc.tol, p.nc.imax, p.sw.newt, p.fuha.lss
%
% See also au2u, fsol, nloopext, getder, resi, stanparam
try Xcont=p.sw.Xcont; catch Xcont=0; end; if Xcont>0; p.Xold=p.X; end % save for bisec 
if (p.fsol.fsol==1 || p.fsol.fsol==3) 
    [u,res,iter,Gu,Glam]=fsol(p,u1); return; end
if p.sw.newt==2; % stopping crit based on stepsize in u, not residual 
    [u,res,iter,Gu,Glam,p]=nloopnleq1(p,u1); return; end 
try almin=p.nc.almin; catch almin=0.2; end % minimal damping, try-catch for backward-comp. 
u=u1; alpha=1; iter=0; 
r=resi(p,u); res0=norm(r,p.sw.norm); % (starting) residual 
Gu=getGu(p,u,r); % derivatives (needed for next tangent, thus always computed) 
skip=0; 
if res0<p.nc.tol % starting residual already small, compute derivatives and return 
   Glam=getGlam(p,u,r); res=res0; skip=1; 
end 
if skip==0; % now start the actual loop
stepok=1; % stepok=1 indicates that step of size alpha decreased residual! 
 % u1(p.nu+p.nc.ilam(1)), pause 
while(res0>p.nc.tol && iter<p.nc.imax && stepok)  
  [upd,p]=p.fuha.lss(Gu,r,p); stepok=0; % iter, pause
  while(stepok==0 && alpha>=almin)       
    au1=u2au(p,u,1) ... % au1=[upde,actaux,lam] 
        -alpha*[upd;0]; % the newton step, no change in primary param. 
    u1=au2u(p,au1,1); 
    r=resi(p,u1); res=norm(r,p.sw.norm);  % residual         
    if res<res0 % good step 
       if Xcont>1; % update p.X, return u in p.up for plotting, set u to zero  
    p.up=u1; [p,u1]=updX(p,u1);  r=resi(p,u1); [Gu,Glam]=getder(p,u1,r); 
       end 
       if(res<res0/4 && alpha<1) % very good step, possibly increase alpha 
           alpha=2*alpha; 
       end 
       stepok=1; u=u1; res0=res;         
       if(p.sw.newt==0) % full Newton, get new derivatives 
          Gu=getGu(p,u,r); % next Jacobian        
       end       
    else alpha=alpha/2; % bad step, try smaller al
    end % if res<res0
  end % while stepok==0 
  iter=iter+1; %alpha=1; 
end   % while res0>p.nc.tol && stepok 
end
% some postprocessing
if(p.sw.newt==0); Glam=getGlam(p,u,r); % Newton, only update Glam 
else [Gu,Glam]=getder(p,u,r); end % chord, update Gu,Glam 
if(p.sw.verb>1); if(alpha<1); % inform user about damping ...
        fprintf('nloop: damp alpha=%g, res=%g\n', alpha,res0); end; end; 
res=res0; % return res 