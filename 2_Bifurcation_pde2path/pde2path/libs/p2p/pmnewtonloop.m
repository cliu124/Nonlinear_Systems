function [uold,tauold,inegp,resp,dsp,pp]=pmnewtonloop(p,p0)
% PMNEWTONLOOP: parallel newtonloops called from pmcont 
%
%  [uold,tauold,inegp,resp,dsp,pp]=pmnewtonloop(p,p0)
m=p.pm.mst; lamd=p.tau(p.nu+p.nc.nq+1);  % lamdot
u0=p.u; ds0=p.sol.ds; tau0=p.tau; p.nc.imax=p.pm.imax;
uold=u0; tauold=tau0; % for bifdec    
iterp=zeros(1,m); resp=zeros(1,m); % for parallel newton loop    
pp=[]; % saves the results like Gu, Glam, lam and u
inegp=-ones(m,length(p.nc.eigref)); % saves eigenvalues  
dsp=zeros(1,m); for i=1:m; dsp(i)=i*ds0; end
try runpar=p.sw.runpar; catch runpar=1; end 
if runpar==1; 
 parfor i=1:m % loop over m predictors 
  au1=u2au(p,p.u,1)+dsp(i)*tau0; u1p=au2u(p,au1,1); % predictor
  moreit=1; res0=100;
  while moreit==1 % for each i, iterate calls to nloop and monitor residual behavior 
    if(p.sw.para==0 || (p.sw.para==1 && abs(lamd)>p.nc.lamdtol)) % fixed lam corrector
        [u1p,resp(i),iter1p,Gup,Glamp]=nloop(p,u1p);meth='nat';
    else [u1p,resp(i),iter1p,Gup,Glamp]=nloopext(p,u1p,dsp(i)); meth='arc'; 
    end
    iterp(i)=iterp(i)+iter1p;
    if resp(i)>p.nc.tol     % no conv yet
      if resp(i)<p.pm.resfac*res0; res0=resp(i);  
      else moreit=0;  end    % stop iteration, no good point
    end
    if resp(i)<p.nc.tol; moreit=0; % convergence to good point
      if(p.sw.spcalc>0);  inegp(i,:)=vspcalc(Gup,p); end % calculate EVals   
      pp(i).Gu=Gup; pp(i).Glam=Glamp; pp(i).u=u1p; 
      pp(i).nu=p.nu; pp(i).nc.ilam=p.nc.ilam;
      pp(i).iter=iterp(i); pp(i).meth=meth; pp(i).sol.ineg=inegp(i,:); 
      if(p.sw.errcheck>0); 
          auxp=p; auxp.u=u1p; err=errcheck(auxp);pp(i).sol.err=err; %HU 
      else pp(i).sol.err=0; 
      end 
    end % if resp(i)<p.tol   
  end % while moreit
 end
else
 for i=1:m % loop over m predictors 
  au1=u2au(p,p.u,1)+dsp(i)*tau0; u1p=au2u(p,au1,1); % predictor
  moreit=1; res0=100;
  while moreit==1 % for each i, iterate calls to nloop and monitor residual behavior 
    if(p.sw.para==0 || (p.sw.para==1 && abs(lamd)>p.nc.lamdtol)) % fixed lam corrector
        [u1p,resp(i),iter1p,Gup,Glamp]=nloop(p,u1p);meth='nat';
    else [u1p,resp(i),iter1p,Gup,Glamp]=nloopext(p,u1p,dsp(i)); meth='arc'; 
    end
    iterp(i)=iterp(i)+iter1p;
    if resp(i)>p.nc.tol     % no conv yet
      if resp(i)<p.pm.resfac*res0; res0=resp(i);  
      else moreit=0;  end    % stop iteration, no good point
    end
    if resp(i)<p.nc.tol; moreit=0; % convergence to good point
      if(p.sw.spcalc>0);  inegp(i,:)=vspcalc(Gup,p); end % calculate EVals   
      pp(i).Gu=Gup; pp(i).Glam=Glamp; pp(i).u=u1p; 
      pp(i).nu=p.nu; pp(i).nc.ilam=p.nc.ilam;
      pp(i).iter=iterp(i); pp(i).meth=meth; pp(i).sol.ineg=inegp(i,:); 
      if(p.sw.errcheck>0); 
          auxp=p; auxp.u=u1p; err=errcheck(auxp);pp(i).sol.err=err; %HU 
      else pp(i).sol.err=0; 
      end 
    end % if resp(i)<p.tol   
  end % while moreit
 end
end