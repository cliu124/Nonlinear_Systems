function [p,iok]=inistep(p)
% INISTEP: return p-structure and ok flag for (attempt of) initial continuation step 
% done by incrementing Parameter and Newton loop.
% This yields first 2 points on branch and first tangent (secant) vector.
%
%  [p,iok]=inistep(p)
%
% Here iok convergence error flag.
%
% See also cont, nloop, stanparam.
fprintf('inistep\n'); try Xcont=p.sw.Xcont; catch; Xcont=0; end; 
iok1=0;iok2=0; [u0,res,iter,~,~,p]=nloop(p,p.u); ineg=-1; muv=[]; %p.sol.ptype=-1;
p.sol.meth='nat';  
if(res<p.nc.tol) 
    if Xcont>0; [p,u0]=updX(p,u0); end % update p.X, return u in p.up, set u to zero    
    p.u=u0; p.sol.res=res; p.sol.iter=iter; 
    if (p.sw.spcalc>0 || p.sw.bifcheck==2) % spectral calc 
        r=resi(p,p.u); Gu=getGu(p,p.u,r);  [ineg,muv]=vspcalc(Gu,p);
    end
    if(p.sw.errcheck>0); p.sol.err=errcheck(p);end 
    if(~isfield(p.sol,'ptype') || p.sol.ptype~=-2) p.sol.ptype=-1; end % initial point, unless from swibra     
    brout=[bradat(p); p.fuha.outfu(p,p.u)]; % userfu to append to bif-branches  
    brplot=brout(length(bradat(p))+p.plot.bpcmp);
    p.branch=extbra(p.branch,brout); % put on branch 
    figure(p.plot.brfig); hold on; % plot point     
    if(any(ineg<=0)) plot(getlam(p),brplot,'*'); else plot(getlam(p),brplot,'+'); end 
    p.sol.ineg=ineg; p.sol.muv=muv;   
    if p.file.smod~=0; p.fuha.savefu(p);  end % save to file 
    dss=0; [p,cstop]=p.fuha.ufu(p,brout,dss);     
    iok1=1;p.file.count=p.file.count+1; 
else fprintf('   - no convergence in zeroth step, lam=%g, res=%g\n',getlam(p),res);
end
u0=p.u; lam0=getlam(p); p.u(p.nu+p.nc.ilam(1))=lam0+p.sol.ds/10; 
%fprintf('   - now lam=%g\n',getlam(p)); 
[u,res,iter,~,~,p]=nloop(p,p.u); 
if(res<p.nc.tol) 
    if Xcont>0; [p,u]=updX(p,u); end % update p.X, return u in p.up, set u to zero
    %fprintf('%4i got second point with lam=%g, res=%g \n',p.file.count,getlam(p),res);
    p.u=u; p.sol.res=res; p.sol.iter=iter; 
    if(p.sw.spcalc>0 || p.sw.bifcheck==2); r=resi(p,u); Gu=getGu(p,u,r); 
        [ineg,muv]=vspcalc(Gu,p); 
    end
    if(p.sw.errcheck>0); p.sol.err=errcheck(p);end 
    brout=[bradat(p); p.fuha.outfu(p,p.u)]; % userfu to append to bif-branches  
    brplot=brout(length(bradat(p))+p.plot.bpcmp);
    p.branch=extbra(p.branch,brout); % put on branch 
    p.sol.ptype=0; % normal point
    figure(p.plot.brfig); hold on; % plot point 
    if(any(ineg<=0)) plot(getlam(p),brplot,'*'); else plot(getlam(p),brplot,'+'); end 
    p.sol.ineg=ineg; p.sol.muv=muv; 
   % if(mod(p.file.count-1,p.file.smod)==0) p.fuha.savefu(p); end  % save to file 
    tau=(u2au(p,p.u,1)-u2au(p,u0,1))/(getlam(p)-lam0);
    p.tau=tau/xinorm(tau,p.sol.xi,p.nc.nq,p.sol.xiq); 
    dss=p.sol.ds/10; [p,cstop]=p.fuha.ufu(p,brout,dss);    
    iok2=1;p.file.count=p.file.count+1;
else fprintf('   - no convergence in first step, lam=%g, res=%g\n',getlam(p),res);
end
iok=iok1*iok2; p.sol.restart=0; 


