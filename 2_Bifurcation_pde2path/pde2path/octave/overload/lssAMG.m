function [x,p]=lssAMG(A,b,p)
% lssAMG: Uses ilupack (www.icm.tu-bs.de/~bolle/ilupack/) for solving Ax=b. 
% MOD for OCTAVE, resp. the ilupack version from https://github.com/fastsolve/ilupack4m
% which mexes to ILUinit etc instead of AMGinit etc 
% 
% [x,p]=lssAMG(A,b,p)
% see also: lss, lssbel, blssbel, belpi, bel
nb=norm(b(:),inf); if nb==0; x=b; return; end   % rhs=0 
if ~isfield(p,'amgopt'); p.amgopt=ILUinit(A,[]); end % very first call of lssAMG 
p.amgopt=ilupparam(p,p.amgopt); % (re)set parameters for AMG 
newprec=0; s=0;   % precon is not yet calculated in this call
if isfield(p.mat,'prec')==1 % check size of preconditioner
   try s=p.mat.prec(1).n; end
end
%s, size(A,1)
if (~isfield(p.mat,'prec') || size(A,1)~=s) % if prec does not ex. or has wrong size 
  if p.sw.verb>1;  disp('computing prec'); end
  if isfield(p.mat,'prec'); % clean memory % problem with prec-size may still give seg-fault
    p.mat.prec=ILUdelete(p.mat.prec); p.mat=rmfield(p.mat,'prec'); 
  end
  t1=tic; [p.mat.prec,p.amgopt]=ILUfactor(A,p.amgopt); t1=toc(t1); newprec=1; % compute prec.
  if p.sw.verb>1; fprintf(' in %g seconds\n', t1); end 
  s=p.mat.prec(1).n; 
end
t1=tic; [x,p.amgopt]=ILUsolver(A,p.mat.prec,p.amgopt,b); t1=toc(t1); 
used_it=max(p.amgopt.niter); 
if p.sw.verb>2; fprintf('ILUsolver in %g seconds, niter=%i \n', t1,used_it); end 
while (used_it>p.amgopt.maxit-1 && newprec==0) % no convergence 
    if newprec==0 % if preconditioner was not updated, do now and try again 
      if p.sw.verb>1; fprintf('\n'); disp('ILUsolver did not converge...new PREC'); end
      p.mat.prec=ILUdelete(p.mat.prec); p.mat=rmfield(p.mat,'prec'); % clean mem 
      p.amgopt=ILUinit(A,[]); p.amgopt=ilupparam(p,p.amgopt); 
      t1=tic; [p.mat.prec,p.amgopt]=ILUfactor(A,p.amgopt); t1=toc(t1); newprec=1; 
      if p.sw.verb>1; fprintf(' in %g seconds\n', t1); end 
      t1=tic; [x,p.amgopt]=ILUsolver(A,p.mat.prec,p.amgopt,b); t1=toc(t1); 
      used_it=max(p.amgopt.niter); 
      if p.sw.verb>2; fprintf('AMGsolver in %g seconds, niter=%i \n', t1, used_it); end       
    end
    if used_it>p.amgopt.maxit-1 % if r is still too large or prec already new 
       choi=0; try choi=p.ilup.noconvhandling; catch; end % default 
       if p.sw.inter>0 
       disp('AMGsolver did not converge. Your options:  '); 
       disp(' 0: stop; 1: decrease droptol, 2: increase maxiter'); 
       choi=asknu('your choice: ',choi); 
       p.ilup.noconvhandling=choi; 
       end
       switch choi; 
           case 0;  error(['AMGsolver did not converge. Stopping.']);  
           case 1;  % decrease droptol and compute new prec 
          if p.ilup.droptol>10*(p.ilup.droptolmin-1e-10)
            p.ilup.droptol=p.ilup.droptol/10; p.ilup.droptolS=p.ilup.droptolS/10; 
            p.amgopt=ilupparam(p,p.amgopt); 
            disp(['decreasing droptol to ' mat2str(p.amgopt.droptol)]); 
            newprec=0; 
          else % clean mem and stop
            p.mat.prec=ILUdelete(p.mat.prec); p.mat=rmfield(p.mat,'prec'); % clean mem  
            error(['AMGsolver did no converge with droptol=' mat2str(p.amgopt.droptol) '. Stopping.']);  
          end
           case 2; % increase maxit and try to solve again 
          if p.amgopt.maxit<(p.ilup.maxitmax+1)/2; 
            p.ilup.maxit=p.ilup.maxit*2; p.amgopt=ilupparam(p,p.amgopt); 
            disp(['solving with maxit=' mat2str(p.amgopt.maxit)]); 
            t1=tic; [x,p.amgopt]=ILUsolver(A,p.mat.prec,p.amgopt,b); t1=toc(t1); 
            fprintf('AMGsolver in %g seconds, niter=%i \n', t1,used_it);
          else % clean mem and stop 
            p.mat.prec=ILUdelete(p.mat.prec); %p.mat=rmfield(p.mat,'prec'); % clean mem  
            error(['AMGsolver did no converge with maxit=' mat2str(p.amgopt.maxit) '. Stopping.']); 
          end
       end
    end
end