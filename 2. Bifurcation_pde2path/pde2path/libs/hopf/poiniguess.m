function p=poiniguess(p,t,u,varargin)
% poiniguess: do iniguess for periodic orbit (indep. of branch-switching)
% set-up data-structures to subseq.call honloop(ext) or (ho)cont
%  
% p=poiniguess(p,t,u,varargin)
%
% t=time-discret, tl points, u=matrix nu x tl
% aux may include
%  nqh=#of Hopf constraints, 
%  ds=stepsize
%  utau=user tangent (if not given, then tau=u is assumed (coming from 0))
%
% spatial discret., FEM matrices etc already in p; essentially mimics hoswibra
aux=[]; if nargin>3; aux=varargin{1}; end 
try nqh=aux.nqh; catch nqh=0; end 
try ds=aux.ds; catch; ds=0.1; end 
try utau=aux.tau; catch; utau=[]; end % user tangent tau 
p=hostanparam(p,aux); % set standard values 
tl=length(t); T=t(end); p.hopf.tl=tl; p.hopf.t=t./T; % rescaled discr.
p.hopf.T=T; nu=p.nu; p.hopf.xi=p.hopf.xif/(nu*tl); 
lam=getlam(p); dlam=0; p.hopf.lam=lam; 
p=resetc(p); brout=[bradat(p); p.fuha.outfu(p,p.u)]; brout(1,1)=0; 
p.hopf.y=u; 
p.hopf.nqh=nqh; 
hoplot(p,p.plot.pfig,p.plot.pcmp,p.hopf.aux); 
figure(p.plot.pfig); title('initial guess');
p.sw.para=4; p.sw.spcalc=0; p=hoMini(p); p.sol.ds=ds; par=p.u(p.nu+1:end); 
p.sol.restart=0; % don't treat as a restart (would need a bit more 'tuning'), 
tau=zeros(1,nu*tl+2+nqh);     % the Hopf tangent vector
if isempty(utau); % fill u-part of tangent with u itself (assuming to be coming from 0)) 
  for i=1:tl; si=(i-1)*nu+1; tau(si:si+nu-1)=u(:,i)/ds;  end;  
else 
  for i=1:tl; si=(i-1)*nu+1; tau(si:si+nu-1)=utau(:,i)/ds;  end;   
end 
tau=tau/honorm(p,tau); p.hopf.tau=tau; 
p.hopf.y0d=sety0dot(p,p.hopf.y,par,p.hopf.T); % phase cond 
%p.branch=brout; p.file.count=0; p.sol.ptype=3; p.fuha.savefu(p);p.file.count=1; 
