function [p,bif]=hobifdetec(p,y,tau,jac,muv1,muv2,ind)
% hobifdetec: Hopf version of bifdetec, para=4 required
try; if ind==p.hopf.oldind; bif=0; return; end; catch; bif=0; return; end; 
try; freeT=p.hopf.freeT; catch; freeT=1; end  % check if T is free (default) 
try; pcsw=p.hopf.pc; catch; pcsw=1; end 
try; nqh2=length(p.hopf.ilam); catch; nqh2=0; end 
try; a=p.u(p.nu+p.hopf.ilam); catch; a=0; end; a1=a; 
if ~p.hopf.indini || p.branch(2,end-1)==1; bif=0; return; end
mucand=0.9; try; mucand=muv1(end-1); end % candidate
if ind>p.hopf.oldind; try; mucand=muv2(1); end; end % if new unstable multiplier(s)
fprintf('possible Bif from Hopf, old-ind=%i, ind=%i, abs(mu)=%g, mu=%g+%gi\n', ...
         p.hopf.oldind, ind, abs(mucand), real(mucand), imag(mucand)); 
bif=1; %q=p; % to later restore
tl=p.hopf.tl; nqh=p.hopf.nqh; na=p.nu+p.nc.nq; ds=p.sol.ds; % just some shorthands 
try; hobisec=p.hopf.bisec; catch; hobisec=3; end 
for i=1:hobisec
    ds=ds/2;  % half stepsize
    fprintf('ds=%g,  ',ds); lam0=p.hopf.lam;  y1=p.hopf.y+ds*v2tom(p,p.hopf.tau);   
    if freeT;  T1=p.hopf.T+ds*p.hopf.tau(na*tl+1); 
        lam1=lam0+ds*p.hopf.tau(na*tl+2);      
        if nqh>0; a0=p.u(p.nu+p.hopf.ilam);   % Hopf with constraints 
          a1=a0+ds*p.hopf.tau(na*tl+3:end)'; % aux-vars predictor, 
        end         
    else % fixed T, nqh>0 also dealt with here 
        T1=p.hopf.T; lam1=lam0+ds*p.hopf.tau(na*tl+1);
        if pcsw
         a0=p.u(p.nu+p.hopf.ilam); 
         a1=a0+ds*p.hopf.tau(na*tl+2:end)'; % aux-vars predictor, aux-tangent at end            
        end
    end     
    [y,T,lam,a,res,iter,A,p]=honloopext(p,y1,T1,lam1,a1,ds);     
    p.hopf.T=T; p.hopf.lam=lam; p.hopf.y=y; p=setlam(p,lam); p.sol.res=res;
    if nqh2>0; p.u(p.nu+p.hopf.ilam)=a; end     
    jac=A(1:p.nu*tl, 1:p.nu*tl); [muv1,muv2,ind1]=floq(p,jac); % new floq
    muc=getmuc(muv1,muv2); fprintf('ga_crit=%g\n',muc); 
    if ind1==p.hopf.oldind; % before BP, do nothing 
    else p.hopf.oldind=ind1; ds=-ds; % behind BP, flip direction 
    end 
    %ds, getlam(p), pause 
    taurhs=[zeros(na*tl+pcsw,1); 1]; % new tangent 
    if nqh>0; taurhs=[taurhs; zeros(nqh,1)]; end 
    [tau,p]=p.fuha.blss(A,taurhs,p); 
    hn=honorm(p,tau); tau=tau/hn; p.hopf.tau=tau';  % new tangent
end % bisections done, now postprocess
p.hopf.indini=0; p.hopf.oldind=ind;
p.hopf.muv1=muv1; p.hopf.muv2=muv2; p.sol.ptype=1; 
if isfield(p.hopf,'auxp'); p.hopf.auxp(p); end % auxiliary plot 
brout=[bradat(p); p.fuha.outfu(p,p.u)];     % userfu to append to bif-branches  
brplot=brout(length(bradat(p))+p.plot.bpcmp);    % y-axis value in bif-figure
p.branch=[p.branch(:,1:end-1) brout]; % put on branch and save 
p.fuha.savefu(p); %p.file.count=p.file.count+1;  
fprintf(['\n BP saved to ' p.file.dir '/bpt' mat2str(p.file.bcount,3) '\n']); 
p.file.bcount=p.file.bcount+1; 
hoplot(p,p.plot.pfig,p.plot.pcmp); figure(p.plot.brfig); hold on; 
plot(getlam(p),real(brplot),'o'); drawnow; 
% older stuff, restore point behind bif; throwing it away seems cleaner 
%p=q; p.branch=[bs, q.branch(:,end)]; p.branch(1,end)=p.branch(1,end)+1; 
%p.file.bcount=p.file.bcount+1; p.file.count=p.file.count+1;  