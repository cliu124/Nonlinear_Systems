function p=ulamcheckho(p)
% ULAMCHECKHO: called by ulamcheck, if hopf setting 
%
%    p=ulamcheck(p)
try; freeT=p.hopf.freeT; catch; freeT=1; end % check if T is free (default) 
targets=[p.usrlam p.nc.lammin p.nc.lammax]; ind=[];
if(size(p.branch,2)>1)
    for j=1:length(targets)
        if(sign(p.branch(4,end)-targets(j))~=sign(p.branch(4,end-1)-targets(j)))
            ind=j; break
        end
    end
end
if(~isempty(ind) && length(p.hopf.tau)>1)
fprintf('computing sol at user defined lambda=%g\n',targets(j)); 
slss=p.fuha.lss; sblss=p.fuha.blss; try; sinnerlss=p.fuha.innerlss; catch; end % save lssolvers, 
%p.fuha.lss=@lss; p.fuha.blss=@lss; p.fuha.innerlss=@lss; % cause otherwise AMG seg-faults
lamd=p.hopf.tau(end); ds=(targets(ind)-p.branch(4,end))/lamd;
y1=p.hopf.y; +ds*v2tom(p,p.hopf.tau);  lam1=p.hopf.lam+ds*lamd; % predictor
if freeT; T1=p.hopf.T+ds*p.hopf.tau(p.nu*p.hopf.tl+1); else; T1=p.hopf.T; end 
[y,T,lam,res,iter,A,q]=honloop(p,y1,T1,lam1);  % corrector
if(res<=p.nc.tol) %if convergence in Newton
   p.file.count=p.file.count+1; % p.branch(1,end)=p.branch(1,end)+1; 
   q.hopf.y=y; q.hopf.lam=lam; q.hopf.T=T; p.branch(1,end)=p.branch(1,end)+1; 
   q.sol.ptype=-3; % flag for usrlam!    
   if(mod(p.file.count,p.file.smod)==0) % if last cont-pt was saved,     
     p.fuha.savefu(p); end  % resave it with label increased by 1         
     ind=-1; 
     if p.sw.para>3 && p.hopf.flcheck>0    
      jac=A(1:q.nu*q.hopf.tl, 1:q.nu*q.hopf.tl); 
      if p.hopf.flcheck==1; [muv1, muv2, ind]=floq(p,jac); 
      else [muv1, muv2, ind]=floqps(p,jac);
      end
      %fprintf('ind=%i', ind); 
      if q.hopf.indini && ind~=q.hopf.oldind
         try mucand=muv1(end-1); catch mucand=0; end % candidate
         if ind>q.hopf.oldind; mucand=muv2(1); end % if new unstable multiplier(s)
         fprintf('possible Bif from Hopf, old-ind=%i, ind=%i, abs(mu)=%g, mu=%g+%gi\n', ...
             q.hopf.oldind, ind, abs(mucand), real(mucand), imag(mucand)); 
      end
      q.hopf.indini=1; q.hopf.oldind=ind;
    end            
   brout=[bradat(q); p.fuha.outfu(q,q.u)];  
   brplot=brout(length(bradat(q))+p.plot.bpcmp); %y-axis value in bif-figure        
   fprintf('%4i %s %s %5.2e %4i  %s %5.2e ', ...
       q.file.count,printcon(getlam(q)),printcon(brplot),res,iter,printcon(ds),T);
   npr=length(q.sw.bprint);
   for i=1:npr; fprintf('%s ',printcon(brout(length(bradat(q))+p.sw.bprint(i)))); end;
   fprintf('\n');
   % put on branch (the new point ends up between the last and the one but last)
   p.branch=[p.branch(:,1:end-1) brout p.branch(:,end)];  
   %p.branch(1,end)=p.branch(1,end)+1; 
   q.branch=p.branch; q.file.count=p.file.count-1;
   p.fuha.savefu(q); % save the solution at the desired lambda to file with label p.file.count
   figure(p.plot.brfig); hold on; try; plot(getlam(q),real(brplot),'sq'); catch; end   % plot point  
else fprintf('no convergence in holoop for lambda = %g \n', targets(ind)); %do nothing (for now)
end
p.fuha.lss=slss; p.fuha.blss=sblss; try; p.fuha.innerlss=sinnerlss; catch; end
end


