function p=cont(p,varargin)
% CONT: main continuation routine 
%
%  p=cont(p)   :  do p.nc.nsteps continuation steps
%  p=cont(p,n) :  do n continuation steps 
%
% (initial) stepsize p.sol.ds, starting from p.u, direction tau, if it's there
q=p; % store initial p in case of failure.
try Xcont=p.sw.Xcont; catch; Xcont=0; end 
if p.sw.sfem~=0 
    if ((p.sw.spcalc~=0 && isempty(p.mat.M)))  
    fprintf('\nNo mass matrix  -- calling p=setfemops(p)\n');
    p=setfemops(p); 
    end
end
if(p.file.pnamesw==1) % set filenames to current input variable name
    [p,ok]=setfn(p,inputname(1));  if ok~=1; return; end
end
% axes-labels
auxdictl=0; if(isfield(p.plot,'auxdict')) auxdictl=length(p.plot.auxdict); end
if(p.nc.ilam(1)<=auxdictl) xax=p.plot.auxdict{p.nc.ilam(1)};  
else xax=['parameter=' mat2str(p.nc.ilam(1))]; end
figure(p.plot.brfig); xlabel(char(xax)); 
yax=mat2str(p.plot.bpcmp);
if(p.plot.bpcmp>0)
  if(p.plot.bpcmp<=auxdictl) yax=p.plot.auxdict{p.plot.bpcmp};  
  else yax=['user branch comp. ' yax]; end
else
  switch p.plot.bpcmp
    case 0; yax='L2-norm';
    case -1; yax='err';
    otherwise; yax=length(bradat(p))+p.plot.bpcmp; yax=['branch comp. ' yax];
  end
end
ylabel(char(yax));

p.time.newton=0; p.time.totst=0; p.time.spec=0; totime=tic; % reset timing 
msteps=p.nc.nsteps;if nargin>1; msteps=varargin{1}; end
if p.sw.para>2; p=hocont(p,msteps); return; end
is=0; % stepcounter 
ineg0=p.sol.ineg; % to save last #neg EVals
switch p.sw.spcont; 
 case 1; fprintf('BP continuation. Use p=bpcontexit(p) to return to normal continuation.\n'); 
 case 2; fprintf('FP continuation. Use p=spcontexit(p) to return to normal continuation.\n'); 
 case 3; fprintf('HP continuation. Use p=hpcontexit(p) to return to normal continuation.\n');
end
p.fuha.headfu(p);   % output user defined header
while is<msteps;    % ********************** continuation loop ****************************
   if (length(p.tau)~=p.nu+p.nc.nq+1)||(p.sol.restart>0) % initial step 
      if(p.file.count>0) p.file.count=p.file.count-1; end
      [p,iok]=inistep(p); if iok==0; p=q; return; end 
      ineg0=p.sol.ineg; 
      plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle); is=is+1; 
      if is>=msteps; p.file.count=p.file.count-1; p.fuha.savefu(p); return; end
   end 

   stime=tic; iter=0; res=10*p.nc.tol; stepok=0; % stepok=1 after successful step 
   dss=p.sol.ds; % save current ds for restart at current lam after mesh-ref-in cfail 
   while stepok==0     % loop to find next point
      lamd=p.tau(p.nu+p.nc.nq+1);
      au1=u2au(p,p.u,1)+p.sol.ds*p.tau; u1=au2u(p,au1,1); % predictor
      %try p.X=predX(p,u1); catch; end 
      ntime=tic; 
      if(p.sw.para==0 || (p.sw.para==1 && abs(lamd)>p.nc.lamdtol)) % fixed lam corrector
          [u1,res,iter,Gu,Glam,p]=nloop(p,u1); p.sol.meth='nat'; 
      else [u1,res,iter,Gu,Glam,p]=nloopext(p,u1,p.sol.ds); p.sol.meth='arc'; % arclength-corrector   
      end
      p.time.newton=p.time.newton+toc(ntime); % newton-loop time, accumulated
      dss=p.sol.ds; % for output below, now stepsize control (re convergence)      
      [p,stepok,u1,res,iter,Gu,Glam]=sscontrol(p,u1,res,iter,Gu,Glam,dss); 
      if(stepok==-1); p=q; return; end; % ABORT cont 
   end              % stepok==0     
   % **********  step accepted    
   p.sol.ptype=0; %  so far normal point, check EVals, bifs, etc 
   if p.sw.spcalc>0 % calculate EVals 
       sptime=tic; ineg0=p.sol.ineg; [p.sol.ineg,p.sol.muv]=vspcalc(Gu,p); 
       p.time.spec=p.time.spec+toc(sptime); % spectral-time, accumulated 
   end
   try secpred=p.sw.secpred; catch; secpred=0; end 
   if secpred % secant instead of tangent 
     ua2=u2au(p,u1,1); ua1=u2au(p,p.u,1); tau1=sign(p.sol.ds)*(ua2-ua1); 
   else  % form extended matrix and compute new tangent
     amat=genamat(p,Gu,Glam,p.tau,p.sol.xi,p.sol.xiq); 
     [tau1,p]=p.fuha.blss(amat,[zeros(p.nu+p.nc.nq,1);1],p); 
   end
   tau1=tau1/xinorm(tau1,p.sol.xi,p.nc.nq,p.sol.xiq);
   fc=p.file.count;    
   if ~isfield(p.sw,'abs'); p.sw.abs=0; end % abs=1 if we are directly after a bifurcation 
   if p.sw.bifcheck>0 && p.sw.abs~=1  % check for bifurcation, unless 1st step after bif
     newds=p.sol.ds; p.sol.ds=dss; % use ds from BEFORE last sscontrol! 
     biftime=tic; [p,bif]=bifdetec(p,u1,tau1,Gu,Glam,ineg0); p.time.bif=p.time.bif+toc(biftime); 
     ineg0=p.sol.ineg;      
     if bif; if p.sw.cdbb==1; u1=p.u; p.sol.ds=newds/2; end % cont directly behind bif                      
     else p.sol.ds=newds; 
     end
   end
   if(p.sw.foldcheck>0 && abs(sign(tau1(p.nu+p.nc.nq+1))-sign(lamd))>1 &&p.sw.abs~=1)  % fold via sgn(lamd)
     newds=p.sol.ds; p.sol.ds=dss; % use ds from BEFORE last sscontrol! 
     p=folddetec(p,u1,tau1); p.sol.ds=newds; 
   end         
   if Xcont>0; [p,u1]=updX(p,u1); % update p.X, return u in p.up, set u to zero
      r=resi(p,u1); [Gu,Glam]=getder(p,u1,r);
   end 
   is=is+(p.file.count-fc);  p.sol.restart=0; p.sw.abs=0;  
   p.u=u1;  p.tau=tau1;    % *********** update u and tau **************
   p.sol.res=res; p.sol.iter=iter; p.sol.lamd=lamd; % store stuff (if e.g. p.fuha.ufu needs iter)
   if(p.sw.errcheck>0); p.sol.err=errcheck(p);  % ERRCHECK and possibly mesh-refinement
     if(p.sol.err>p.nc.errbound && p.nc.errbound>0)                            
       if(p.sw.errcheck==1  || p.sw.bcper~=0) % just give warning! 
         fprintf('   - err.est.=%g>errbound=%g. Consider mesh-refinement.\n', p.sol.err, p.nc.errbound);          
       end 
       if(p.sw.errcheck==2 && p.sw.bcper==0); % adapt mesh 
           fprintf('   - err.est.=%g>errbound=%g. Adapting mesh\n', p.sol.err, p.nc.errbound); 
           p=meshadac(p,'eb',p.nc.errbound); 
       end 
       if(p.sw.errcheck>2) && p.sw.bcper==0; p=meshref(p,'eb',p.nc.errbound); end % refine mesh 
     end %p.sol.err>p.nc.errbound  
   end %p.sw.errcheck>0  
   % Check for various actions
   if (mod(p.file.count,p.nc.amod)==0 && p.file.count~=0); % adapt mesh ignoring errbound 
       fprintf('   - adapting mesh\n'); p=meshadac(p,'eb',0);  
   end    
   brout=[bradat(p); p.fuha.outfu(p,p.u)];          % userfu to append to bif-branches  
   brplot=brout(length(bradat(p))+p.plot.bpcmp);    %y-axis value in bif-figure
   p.branch=extbra(p.branch,brout); %[p.branch brout];   % put on branch 
   figure(p.plot.brfig); grid on; hold on;                   % plot point with according symbol
   if p.sol.ineg<=0; plot(getlam(p),real(brplot),'*b'); drawnow; 
   else plot(getlam(p),real(brplot),'+b'); drawnow; end 
   p.time.totst=p.time.totst+toc(stime);            % total step time (accumulated) 
   [p,cstop]=p.fuha.ufu(p,brout,dss);                % user function, typically printout    
  if(p.file.count>0 && mod(p.file.count,p.file.smod)==0) % save to file 
      p.fuha.savefu(p); end
   if(mod(p.file.count,p.plot.pmod)==0); 
       plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle);
   end % plot sol       
   p.file.count=p.file.count+1; is=is+1; 
   if p.file.count>p.nc.ntot; cstop=1; end; 
   if isfield(p.mat,'prec')
       try
   if(mod(p.file.count-1,p.ilup.ilun)==0); % force update of prec
     disp('new prec');
     p.mat.prec=AMGdelete(p.mat.prec); p.mat=rmfield(p.mat,'prec');   
     [p.mat.prec,p.amgopt]=AMGfactor(Gu,p.amgopt);
   end
       end
   end
   if cstop==1; break; end                          % p.fuha.ufu returned stop! 
end % while is<msteps
p.file.count=p.file.count-1; % some postprocessing, i.e., save the last point 
if(mod(p.file.count,p.file.smod)~=0 && p.file.smod~=0); p.fuha.savefu(p); end % save last point with adjusted counter 
p.file.count=p.file.count+1; p.time.tot=toc(totime); 
if(p.time.timesw>0); fprintf('Timing: total=%g, av.step=%g, av.Newton=%g, av.spcalc=%g\n',...
      p.time.tot,p.time.totst/is,p.time.newton/is,p.time.spec/is);
end
