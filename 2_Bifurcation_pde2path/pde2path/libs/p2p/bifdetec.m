function [p,bif]=bifdetec(p,u1,tau1,Gu,Glam,ineg0)
% BIFDETEC: Called from cont to detect bifurcation 
% If bif occurs, locate bif-point by bisection and save it. 
% If p.sw.bifcheck=1 use LU-decomp to calc  sign(det(A)) 
%                 =2 count neg EVals near shifts p.nc.eigref(j) 
% 
%  p=bifdetec(p,u1,tau1,Gu,Glam,ineg0,inegn)
%
%  p.u, p.sol.deta, ineg0    from previous continuation point data  
%  u1, lam1, tau1, Gu, Glam  at new point.
%
global p2pglob; 
bif=0; deta0=p.sol.deta; ds=p.sol.ds; xi=0.5; sds=sign(ds); 
Gubb=Gu; % use balanced xi for det(A) 
amat=genamat(p,Gu,Glam,tau1,xi,xi); muv=0; 
try Xcont=p.sw.Xcont; catch; Xcont=0; end; 
switch p.sw.bifcheck
    case 1; try; [L,U,P,Q]=lu(amat); j0l=1;  p.sol.deta=full(prod(sign(diag(U)))*det(P)*det(Q));
        catch; [L,U]=lu(amat); j0l=1;  p.sol.deta=full(prod(sign(diag(U)))); end 
    if (p.sol.restart>0 || abs(p.sol.deta-deta0)<1.1); return; end  % restart, or no bif, do nothing 
    case 2; bif=0;  j0l=[]; 
        [ineg1,muv]=vspcalc(Gu,p); % compute ineg at different om's      
        for j=1:length(p.nc.eigref)
        try;  if (ineg1(j)~=ineg0(j) && abs(real(muv(j)))<p.nc.mu1); 
                 bif=1; j0l=[j0l, j]; end; 
        catch; end 
        end
        if (p.sol.restart>0 || bif==0); p.sol.ineg=ineg1; return; end  % restart, or no bif, do nothing 
end
bif=1;  % Bifurcation detected, possibly cross-check with ineg;
if(p.sw.verb>0)
  if p.sw.bifcheck==1; fprintf('   bifurcation between %g and %g\n',getlam(p),getlam(p,u1));  
  else fprintf('   %i possible bifurcation between %g and %g, om=',length(j0l),getlam(p),getlam(p,u1));  
    for j=1:length(j0l)-1; fprintf('%g,', imag(p.nc.eigref(j))); end; 
    fprintf('%g\n',imag(p.nc.eigref(j0l(end)))); 
  end 
end
% Now localize! 
u0s=p.u; tau0s=p.tau; dsorg=p.sol.ds; 
if Xcont>0; Xsave=p.X; if Xcont==3; p.X=p.Xold; end; end  
ubb=u1;  % ubb (behind bif) used for next cont if sw.cdbb=1 %AM
for j0=j0l % loop over shifts
   u0=u0s; tau0=tau0s; ds=ds/2; bisecc=0; wb=0; % start bisection, wb stands for was behind bif
  lam0=getlam(p,u0); lam1=getlam(p,u1); % lam0, lam1
   while (bisecc<p.nc.bisecmax && abs(ds)>p.nc.dsminbis)        
     au0=u2au(p,u0,1); au1=u2au(p,u1,1); % get active variables
     if p.sw.bifloc==0 % tangent predictor 
       aun=au0+ds*tau0; % select entries of tau0 into full u = [upde;par]
     else if p.sw.bifloc==1; aun=au0+0.5*(au1-au0); % secant predictor 
          else aun=0.25*(3*au0+au1)+0.5*ds*tau0;  end % quadratic predictor
     end
     un=au2u(p,aun,1); u0=au2u(p,au0,1); % merge with auxiliary variables
     if(p.sw.verb>1) %fprintf(' checking lam=%g ...',aun(p.nu+p.nc.nq+1)); end 
     fprintf(' checking lam=%f ...',aun(p.nu+p.nc.nq+1)); end 
     if (p.sw.para==0 || (p.sw.para==1 && abs(tau0(p.nu+p.nc.nq+1))>p.nc.lamdtol))
      [un,res,iter,Gu,Glam,p]=nloop(p,un);  % fixed lam corrector
     else [un,res,iter,Gu,Glam,p]=nloopext(p,un,ds,u0); % arclength corrector
     end
     if(res<p.nc.tol) % If step was (rather) OK, then 
      amat=genamat(p,Gu,Glam,tau0,p.sol.xi,p.sol.xiq);
      [tau1,p]=p.fuha.blss(amat,[zeros(p.nu+p.nc.nq,1);1],p);
      tau1=tau1/xinorm(tau1,p.sol.xi,p.nc.nq,p.sol.xiq); % tau at new point  
      amat=genamat(p,Gu,Glam,tau1,xi,xi); % at new point 
      if p.sw.bifcheck==1; 
         try;  [L,U,P,Q]=lu(amat); detab=full(prod(sign(diag(U)))*det(P)*det(Q));
         catch; [L,U]=lu(amat); detab=full(prod(sign(diag(U)))); end; 
        if (p.sw.verb>1); fprintf(' ok, sign(det(A))=%i\n',detab); end 
        
        if detab==deta0; ds=ds/2; % u_n still before bif, reset u0, keep u1;          
        else u1=u0; deta0=detab; ds=-ds/2; % u_n is behind bif., reset u1 to u0, i.e., flip direction:        
            if Xcont>0;  % update p.X, return u in p.up, set u to zero                 
         [p,un]=updX(p,un); r=resi(p,un); [Gu,Glam]=getder(p,un,r);                     
            end
        end
      else [ineg(j0),muv,V]=spcalc(Gu,p,j0);     
        if (p.sw.verb>1); fprintf(' ok, ineg=%i\n',ineg(j0)); end
        if ineg(j0)==ineg0(j0); ds=sds*abs(ds)/2; % u_n is before bif, sds~sign(ds0)          
            if wb==1; u1=u0; end % set u1 to u0 if step before was behind bif (wb)          
            wb=0; 
        else ds=-sds*abs(ds)/2; % u_n is behind bif, sds=sign(ds0)            
            if wb==0; u1=u0;  end; wb=1; % set u1 to u0 if step before was before bif (wb)
        end
      end
      if sign(ds)~=sign(dsorg); %
          ubb=un; 
        if p.sw.cdbb==1; if p.sw.bifcheck==2; ineg1=ineg; if Xcont>0; Xsave=p.X; end;  else Gubb=Gu; end; end % AM
                                                             % Zwei zeilen tiefer als vorher 
      end % u1 is "behind" (in the org sense) bif 
      tau0=tau1;                  
      if Xcont>0;  % update, but only behind bif point       
         [p,un]=updX(p,un); r=resi(p,un); [Gu,Glam]=getder(p,un,r);                   
      end                            
      u0=un;  % set tau, u0 for the next step in bisec       
    else % if res<2*p.nc.tol
      if(p.sw.verb>1) fprintf('No convergence, localization might be poor ...\n'); end
      ds=3*ds/4; % might have hit sing point, try a smaller step in predictor
    end % if res<2*p.nc.tol 
    bisecc=bisecc+1;
  end % localization done 
  if Xcont>0; % final updX, HU    
      [p,un]=updX(p,un); r=resi(p,un); [Gu,Glam]=getder(p,un,r); 
  end; %p.X=Xsave;
  if p.sw.bifcheck==2 && abs(real(muv(1)))>p.nc.mu2 
     fprintf('mu_r=%g, mu_i=%g, no convergence\n',real(muv(1)),imag(muv(1)));  
     bif=0; try p.sol.ineg(j0)=ineg(j0); catch; end;  % update ineg to avoid double checks
   %  p.sol.ineg(:)=ineg1(:); % 20-01
  else 
      hopf=0; try hptol=p2pglob.hptol; catch; hptol=1e-4; end 
  if p.sw.bifcheck>1
      fprintf(' mu_r=%g, mu_i=%g \n',real(muv(1)),imag(muv(1)));
      if abs(imag(muv(1)))>hptol; hopf=1; end; 
  end
  p.sol.j0=j0; % eigref-nr 
%------------------------------------------------------- postprocessing!
  p.u=u0; p.tau=tau0; % and store current values for save 
  if (p.sw.spcalc>0 || p.sw.bifcheck==2) % ineg at bifpoint
      [p.sol.ineg,p.sol.muv]=vspcalc(Gu,p); end  
  p.sol.ptype=1; if hopf; p.sol.ptype=3; end  % flag type of bifurcation point
  if(p.sw.errcheck>0); p.sol.err=errcheck(p); end 
  if (p.sol.ptype==1 && p.sw.bifcheck==2) % to exclude FPs, check if Glam\in R(Gu)=N(Gu)^T
    if p.sw.eigsstart==1; vs=size(Gu,1); p.sw.evopts.v0=ones(vs,1)/vs; 
    else; try p.sw.evopts=rmfield(p.sw.evopts,'v0'); catch; end; end 
    [phi,mu]=eigs(Gu,1,0,p.sw.evopts); phi=phi/norm(phi);  
    [psi,mu]=eigs(Gu',1,0,p.sw.evopts); psi=psi/(phi'*psi); 
    del=Glam'*psi; fprintf('<phi,psi>=%g,',del); 
    try foldtol=p.nc.foldtol; catch; foldtol=0.01; end 
    if abs(del)>foldtol; fprintf(' Fold\n'); p.sol.ptype=2; 
        if p.sw.foldcheck==0; p.u=ubb; if Xcont~=0; p.X=Xsave; end; return; end ; %
    else; fprintf('BP\n'); 
    end 
  end
  brout=[bradat(p); p.fuha.outfu(p,u0)]; % userfu to append to bif-branches  
  brplot=brout(length(bradat(p))+p.plot.bpcmp);
  p.branch=extbra(p.branch,brout); p.fuha.savefu(p);   % put on branch and save
  switch p.sol.ptype
      case 1; fname=[p.file.bpname,sprintf('%i',p.file.bcount),'.mat']; ptype='BP';
          p.file.bcount=p.file.bcount+1;
          fprintf(['%4i  %7.5e (' ptype ', saved to %s) bisection steps %i, last ds %g\n'],...
               p.file.count,getlam(p),fname,bisecc,ds);
          figure(p.plot.brfig); plot(getlam(p),brplot,'o');
      case 2;
          fname=[p.file.fpname,sprintf('%i',p.file.fcount),'.mat']; ptype='FP';
          p.file.fcount=p.file.fcount+1;
          fprintf(['%4i  %7.5e (' ptype ', saved to %s) bisection steps %i, last ds %g\n'],...
               p.file.count,getlam(p),fname,bisecc,ds);
          figure(p.plot.brfig); plot(getlam(p),brplot,'x');
      case 3; fname=[p.file.hpname,sprintf('%i',p.file.hcount),'.mat']; ptype='HP';
          p.file.hcount=p.file.hcount+1;
          fprintf(['%4i  %7.5e (' ptype ', saved to %s) bisection steps %i, last ds %g\n'],...
               p.file.hcount,getlam(p),fname,bisecc,ds);
          figure(p.plot.brfig); plot(getlam(p),brplot,'d');
  end
  if((p.file.count>0) && (mod(p.file.count,p.file.smod)==0)) % save to file 
    fname=[p.file.pname,sprintf('%i',p.file.count),'.mat']; % save(fname,'p'); 
   p.fuha.savefu(p,0); fprintf('   also saved to %s\n',fname);
  end 
  p.u=ubb;  if Xcont>0; p.X=Xsave; end, %p.X=p.Xold; % return point after bif but close to bifpoint!  %AM
  if p.sw.bifcheck==2; p.sol.ineg=ineg1;
  else if p.sw.spcalc>0; [p.sol.ineg,p.sol.muv]=vspcalc(Gubb,p);end; 
  end
  p.file.count=p.file.count+1; p.sol.ptype=0; 
  end
end  % for j0=j0l (loop over om_j's) 