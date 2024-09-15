function p=findbifm(p,varargin)
% find bifpoint p1 starting at p such that abs(p1.ineg-p.sol.ineg)=2; 
% if called with varargin{1}=:nbp; try to find nbp bifpoints
% At most p.nc.nsteps steps are done altogether
nbp=1; if nargin>1; nbp=varargin{1}; end; dss=p.sol.ds; dsmaxs=p.nc.dsmax; 
% if cont has not yet been run make just one small step
if p.sol.restart==1; p.file.pnamesw=0; p.sol.ds=1e-10; p=cont(p,1); p.sol.ds=dss; end
% saving parameters, which will change locally to store back below
smod=p.file.smod; bif=p.sw.bifcheck; pnamesw=p.file.pnamesw; 
headfu=p.fuha.headfu; ufu=p.fuha.ufu; timesw=p.time.timesw; pmod=p.plot.pmod;
r=resi(p,p.u); Gu=getGu(p,p.u,r); [p.sol.ineg,p.sol.muv]=spcalc(Gu,p); 
ineg0=p.sol.ineg; % number of negative eigenvalues of starting point
fprintf('Starting findbifm, %i negative EVals at lam=%g\n', p.sol.ineg, getlam(p));
n=0; k=0; st=0; %ft=p;
while (k<nbp && st==0) % loop over nbp desired bif-points 
  k=k+1; af=0; p.file.smod=0; p.file.pnamesw=0;
  p.fuha.headfu=@findbifheadfu;  p.fuha.ufu=@findbiffu; % no saving/output
  p.sw.bifcheck=0; % the bifpoint will we calc. at the end 
  p.time.timesw=0; p.plot.pmod=0; p.sol.ds=dss; p.nc.dsmax=dss; 
  while (st==0 && abs(p.sol.ds)>p.nc.dsminbis) % find next bifpoint
    ft=p; ft=cont(ft,1); n=n+1; 
    if(getlam(ft)<p.nc.lammin); st=2; fprintf('  lam<lammin, stopping\n');end
    if(getlam(ft)>p.nc.lammax); st=2; fprintf('  lam>lammax, stopping\n');end
    if(n>p.nc.nsteps); st=2; fprintf('  n>p.nc.nsteps, stopping\n');end
    if st<2
       r=resi(ft,ft.u);Gu=getGu(ft,ft.u,r);
       [ft.sol.ineg,ft.sol.muv]=spcalc(Gu,ft); 
      if(p.sw.verb>0); fprintf('   - lam=%e, ineg=%i\n',getlam(ft),ft.sol.ineg); end    
      if ft.sol.ineg==ineg0 % index didn't change, next step! 
        p=ft; 
        if af==1; p.sol.ds=p.sol.ds/2; end % index-change already found, decrease ds        
      else % index changed 
        if abs(ft.sol.ineg-ineg0)>2; af=1;   % index-change by more than 2 found           
           if abs(p.sol.ds)>p.nc.dsminbis; p.sol.ds=p.sol.ds/2; end; % decrease ds 
        else st=1; % point found! 
        end
      end
    end
  end   % while find next bifpoint
  % now bisection, starting with current p.sol.ds, 
  fprintf('starting bisec\n'); 
  u0=p.u; u1=ft.u; tau0=p.tau; ds=p.sol.ds/2; bisecc=0; % start bisection 
  while (bisecc<p.nc.bisecmax && abs(ds)>p.nc.dsminbis)  
    au0=u2au(p,u0,1); au1=u2au(p,u1,1); % get active variables
    if p.sw.bifloc==0 % tangent predictor 
       aun=au0+ds*tau0; % select entries of tau0 into full u = [upde;par]
    else if p.sw.bifloc==1;  aun=au0+0.5*(au1-au0);  % secant predictor 
        else aun=0.25*(3*au0+au1)+0.5*ds*tau0; end % quadratic predictor
    end
    un=au2u(p,aun,1); u0=au2u(p,au0,1); % merge with auxiliary variables     
    if(p.sw.para==0 || (p.sw.para==1 && abs(tau0(p.nu+p.nc.nq+1))>p.nc.lamdtol)) % fixed lam corrector
      [un,res,iter,Gu,Glam]=nloop(p,un); 
    else [un,res,iter,Gu,Glam]=nloopext(p,un,ds,u0); % arclength corrector
    end
    p.u=un;lam=getlam(p); 
    if(res<2*p.nc.tol) % If step was OK, then 
       r=resi(p,p.u);Gu=getGu(p,p.u,r);
       [ineg,muv]=spcalc(Gu,p); %fprintf('ds=%e, lam=%e, mu=%e\n',ds,lam,muv(1));
       fprintf('ds=%e, lam=%e, ineg=%i\n',ds,lam,ineg);
       if ineg==ineg0;  ds=ds/2; % u_n still before bif, reset u0, keep u1;        
       else u1=u0; ineg0=ineg; ds=-ds/2;  % u_n is behind bif., flip direction:
       end
       u0=un; amat=genamat(p,Gu,Glam,tau0,p.sol.xi,p.sol.xiq);
       tau1=amat\[zeros(p.nu+p.nc.nq,1);1]; 
       tau1=tau1/xinorm(tau1,p.sol.xi,p.nc.nq,p.sol.xiq);
       tau0=tau1; % set tau, u0,u1 for the next step in bisec 
       bisecc=bisecc+1; 
    else % try different stepsize
    end
  end
  p.u=un; p.branch=[p.branch [bradat(p);p.fuha.outfu(p,un)]]; 
  p.file.smod=smod; p.fuha.headfu=headfu; p.fuha.ufu=ufu;
  p.time.timesw=timesw; p.plot.pmod=pmod; p.sol.ptype=1; 
  fname=[p.file.bpname,sprintf('%i',p.file.bcount),'.mat']; 
  p.fuha.savefu(p); fprintf(['saved to ' fname '\n']); 
  brout=[bradat(p); p.fuha.outfu(p,u0)]; % userfu to append to bif-branches  
  brplot=brout(length(bradat(p))+p.plot.bpcmp);
  figure(p.plot.brfig); plot(getlam(p),brplot,'o');
  
  st=0; p.u=ft.u;
  p.file.bcount=p.file.bcount+1; ineg0=ft.sol.ineg; p.sol.ds=dss;
end % loop over nb (desired) bifpoints
% set remaining parameters back to original
p.sw.bifcheck=bif; p.file.pnamesw=pnamesw; %p.sol.ds=dss; 
p.nc.dsmax=dsmaxs; 
end % findbif

% inner functions
function [p,cstop]=findbiffu(p,brout,ds)
    cstop=0;
end

function findbifheadfu(p)
end
