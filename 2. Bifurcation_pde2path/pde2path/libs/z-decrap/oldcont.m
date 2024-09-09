function p=cont(p,varargin)
% CONT: main continuation routine for PDE -div(c grad u)+au-f-b grad u=0, 
%  and possible auxiliary equations. Here c,a,f,b are encoded in p.fuha.G, BC encoded in p.fuha.bc 
%
%  p=cont(p)   :  do p.nc.nsteps continuation steps
%  p=cont(p,n) :  do n continuation steps 
%
% (initial) size p.sol.ds, starting from p.u, direction tau, if it's already there
%
% See also: inistep, nloop, genamat, nloopext, spcalc, sscontrol, bifdetec, 
%          folddetec, meshref, meshadac, stanparam.
q=p; % store initial p in case of failure.
if(length(p.nc.ilam)~=p.nc.nq+1)
    fprintf('\nLength of p.nc.ilam does not match equations plus primary parameter!\n');
return; end
if( (p.sw.spcalc~=0 && isempty(p.mat.M)) || (p.sw.sfem~=0 && (isempty(p.mat.M) || isempty(p.mat.K))) ) 
    fprintf('\nNo mass or stiffness matrix  -- calling p=setfemops(p)\n');
    p=setfemops(p);
end
if(p.file.pnamesw==1) % set filenames to current input variable name
    [p,ok]=setfn(p,inputname(1));
    if ok~=1; return; end
end
xax=mat2str(p.nc.ilam(1)); yax=mat2str(p.plot.bpcmp);
figure(p.plot.brfig); xlabel(char(['primary parameter = aux. var. ' xax])); 
if(p.plot.bpcmp>0) 
    ylabel(char(['user branch comp. ' yax]));
elseif(p.plot.bpcmp==0) 
    ylabel('L2-norm of u');
else
    yax=length(bradat(p))+p.plot.bpcmp; ylabel(char(['branch comp. ' yax]));
end
msteps=p.nc.nsteps;if nargin>1; msteps=varargin{1}; end
p.time.newton=0; p.time.totst=0; p.time.spec=0; totime=tic; % reset timing 
is=0; % stepcounter 
ineg0=-1; % to save last #neg EVals
if(p.sw.spcont==1) 
    fprintf('Branch point continuation. Use p=spcontexit(p) to return to normal continuation.\n'); 
elseif(p.sw.spcont==2) 
    fprintf('Fold point continuation. Use p=spcontexit(p) to return to normal continuation.\n'); 
end
p.fuha.headfu(p);         % output user defined header
while is<msteps      % continuation loop 
   if (length(p.tau)~=p.nu+p.nc.nq+1)||(p.sol.restart>0) % initial step 
      if(p.file.count>0) p.file.count=p.file.count-1; end
      [p,iok]=inistep(p); if iok==0; p=q; return; end 
      plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle); is=is+1; 
      if is>=msteps; p.file.count=p.file.count-1; p.fuha.savefu(p); return; end
   end
   stime=tic; iter=0; res=10*p.nc.tol; stepok=0; % stepok=1 after successful step 
   dss=p.sol.ds; % save current ds for restart at current lam after mesh-ref-in cfail 
   while stepok==0     % loop to find next cont point
      lamd=p.tau(p.nu+p.nc.nq+1);
      au1=u2au(p,p.u,1)+p.sol.ds*p.tau; u1=au2u(p,au1,1); % predictor
      ntime=tic;
      if(p.sw.para==0 || (p.sw.para==1 && abs(lamd)>p.nc.lamdtol)) % fixed lam corrector
          [u1,res,iter,Gu,Glam]=nloop(p,u1); p.sol.meth='nat'; 
      else % arclength-corrector   
          [u1,res,iter,Gu,Glam]=nloopext(p,u1,p.sol.ds); p.sol.meth='arc'; 
      end
      p.time.newton=p.time.newton+toc(ntime); % newton-loop time, accumulated
      ds=p.sol.ds; % for output below, now stepsize control (re convergence)      
      [p,stepok,u1,res,iter,Gu,Glam]=sscontrol(p,u1,res,iter,Gu,Glam,dss); 
      if(stepok==-1); p=q; return; end; % ABORT cont 
   end              % stepok==0 
   p.sol.ptype=0; % so far normal point
   if p.sw.spcalc>0 % calculate EVals 
       sptime=tic; ineg0=p.sol.ineg;
       [p.sol.ineg,p.sol.muv]=spcalc(Gu,p); 
       p.time.spec=p.time.spec+toc(sptime); % spectral-time, accumulated 
   else p.sol.ineg=-1; 
   end
   % form extended matrix and compute new tangent
   amat=genamat(p,Gu,Glam,p.tau,p.sol.xi,p.sol.xiq); 
   tau1=p.fuha.blss(amat,[zeros(p.nu+p.nc.nq,1);1],p); 
   tau1=tau1/xinorm(tau1,p.sol.xi,p.nc.nq,p.sol.xiq); 
   fc=p.file.count;
   if(p.sw.bifcheck>0)     % check for bifurcation via |det(A)|
     newds=p.sol.ds; p.sol.ds=ds; % use ds from BEFORE last sscontrol! 
     biftime=tic; p=bifdetec(p,u1,tau1,Gu,Glam,ineg0); 
     p.time.bif=p.time.bif+toc(biftime); p.sol.ds=newds; 
   end
   if(p.sw.foldcheck>0 && abs(sign(tau1(p.nu+p.nc.nq+1))-sign(lamd))>1)  % fold via sgn(lamd)
     newds=p.sol.ds; p.sol.ds=ds; % use ds from BEFORE last sscontrol! 
     p=folddetec(p,u1,tau1); p.sol.ds=newds; 
   end
   is=is+(p.file.count-fc);
   p.sol.restart=0; p.u=u1; p.tau=tau1;     % step accepted! 
   p.sol.res=res; p.sol.iter=iter; p.sol.lamd=lamd; % store stuff (if e.g. p.fuha.ufu needs iter)
   if(p.sw.errcheck>0); p.sol.err=errcheck(p);   % ERRCHECK and possibly mesh-refinement
     if(p.sol.err>p.nc.errbound && p.nc.errbound>0)                            
       if(p.sw.errcheck==1  || p.sw.bcper~=0) % just give warning! 
         fprintf('   - err.est.=%g>errbound=%g. Consider mesh-refinement.\n', p.sol.err, p.nc.errbound);          
       end 
       if(p.sw.errcheck==2 && p.sw.bcper==0); % adapt mesh 
           fprintf('   - err.est.=%g>errbound=%g. Adapting mesh\n', p.sol.err, p.nc.errbound); 
           p=meshadac(p,'eb',p.nc.errbound); 
       end 
       if(p.sw.errcheck>2) && p.sw.bcper==0; % refine mesh 
          p=meshref(p,'eb',p.nc.errbound); 
       end          
     end %p.sol.err>p.nc.errbound  
   end %p.sw.errcheck>0  
   % Check for various actions
   if (mod(p.file.count,p.nc.amod)==0 && p.file.count~=0); % adapt mesh ignoring errbound 
       fprintf('   - adapting mesh\n'); 
       p=meshadac(p,'maxt',p.mesh.maxt,'eb',0); 
   end 
   brout=[bradat(p); p.fuha.outfu(p,p.u)];          % userfu to append to bif-branches  
   brplot=brout(length(bradat(p))+p.plot.bpcmp);    %y-axis value in bif-figure
   p.branch=[p.branch brout];                       % put on branch 
   if(p.file.count>0 && mod(p.file.count,p.file.smod)==0) % save to file 
       p.fuha.savefu(p); 
   end
   figure(p.plot.brfig); hold on;                   % plot point with according symbol
   if p.sol.ineg<=0; plot(getlam(p),real(brplot),'*'); 
   else plot(getlam(p),real(brplot),'+'); end 
   p.time.totst=p.time.totst+toc(stime);            % total step time (accumulated) 
   [p,cstop]=p.fuha.ufu(p,brout,ds);                % user function, typically printout 
   if(mod(p.file.count,p.plot.pmod)==0); plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle); end % plot sol
   p.file.count=p.file.count+1; is=is+1; 
   if p.file.count>p.nc.ntot; cstop=1; end; 
   if(cstop==1) break; end                          % p.fuha.ufu returned stop! 
end % while is<msteps

% some postprocessing, i.e., save the last point 
p.file.count=p.file.count-1; 
if(mod(p.file.count,p.file.smod)~=0&&p.file.smod~=0); p.fuha.savefu(p); end % save last point with adjusted counter 
p.file.count=p.file.count+1; 
p.time.tot=toc(totime); 
if(p.time.timesw>0) 
  fprintf('Timing: total=%g, av.step=%g, av.newton=%g, av.spcalc=%g\n',...
      p.time.tot,p.time.totst/is,p.time.newton/is,p.time.spec/is);
end
