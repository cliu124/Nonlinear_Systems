function p=hoswibra(dir,fname,ds,para,varargin)
% HOSWIBRA: branch switching from Hopf points
% 
%  p=hoswibra(dir,fname,ds,para,newdir,aux)
%
%  additional arguments aux include: 
%   z(1,2) for double Hopf: coefficients to mix the two EVecs at bif for the predictor 
%   bpcmp: new branch-compo for plotting
%   xi:    weight for Hopf arclength
%   freeT: if 0 then used fixed period T (needs to free another parameter!)
%   
%   for Newton-Krylov-shooting (NKS, para=6): 
%   nt: # discrt.points 
%   iT: index of T in parameter vector (only for NKS) 
if nargin>5; aux=varargin{2}; ndir=varargin{1}; else; aux=0; end 
try freeT=aux.freeT; catch freeT=1; end; % check if fixed T is desired (freeT=0)
p=loadp(dir,fname); nu=p.nu; p.hopf.freeT=freeT; 
if(~isempty(varargin)); [p,ok]=setfn(p,varargin{1}); 
    if ~ok; fprintf('could not change dir-name!\n'); return; end
else fprintf('warning: problem directory unchanged.\n'); end
p=hostanparam(p,aux); % set standard values, also taking aux into account 
nqh=p.hopf.nqh; tl=p.hopf.tl;
p.hopf.dsw=0; try p.hopf.dsw=aux.dsw; end 
if isfield(aux,'bpcmp'); p.plot.bpcmp=aux.bpcmp; end
pcheck=1; if isfield(aux,'pcheck'); pcheck=aux.pcheck; end; 
p.hopf.xi=p.hopf.xif/(nu*tl); % seems most reasonable
if isfield(aux,'xi'); p.hopf.xi=aux.xi; end; 
lam=getlam(p); p=resetc(p); brout=[bradat(p); p.fuha.outfu(p,p.u)]; brout(1,1)=0; 
p.branch=brout; p.file.count=0; p.sol.ptype=3; p.fuha.savefu(p);p.file.count=1; 
r=resi(p,p.u); Gu=getGu(p,p.u,r); [ineg,muv,V]=spcalc(Gu,p,p.sol.j0); 
try; fprintf(['leading EVals:' num2str(muv(1:4)) '\n']); 
catch fprintf(['leading EVals:' num2str(muv(1:end)) '\n']); end 
evnr=1; om=imag(muv(evnr)); mu0=real(muv(evnr)); 
idx1=abs(real(muv))<p.nc.mu2; muv=muv(idx1); V=V(:,idx1); % extract central evals 
[mus,idx]=sort(abs(imag(muv)-om));  % sort those close to om_0 to front 
muv=muv(idx); V=V(:,idx); muv', evnr 
phi=V(:,evnr); 
phi=phi/norm(phi,2); 
if isfield(aux,'dlam') % user guess for normal form coeffis 
    dlam=aux.dlam; al=10; 
    if isfield(aux,'al'); al=aux.al; fprintf('user set alpha=%g\n',al); end; 
    fprintf('user set dlam=%i, al=%g\n',dlam,al); 
else  % compute normal form coeffis 
  hodel=1e-4; if isfield(aux,'hodel'); hodel=aux.hodel; end 
  [dlam,al]=hogetnf(p,Gu,om,mu0,phi,hodel);
end 
fprintf('real(mu)=%g, imag(mu)=%g\n',mu0,imag(muv(1)));  
fprintf('lam=%g, om=%g,  NF-coeffis: dlam=%i, alpha=%g\n',lam, om, dlam,al); 
if isfield(aux,'nqnew'); p.nc.nq=aux.nqnew; end % reset nq for Mu_t=-resi, 
tl=p.hopf.tl; t=linspace(0,1,tl); 
p.hopf.t=t; p.hopf.T=2*pi/abs(om); p.hopf.y=p.u(1:nu)*ones(1,tl); 
if ~isfield(aux,'z') % standard initial guess
  ig=ds*al*real(phi(1:nu)*exp(-2i*sign(om)*pi*t));  
else % initial guess, version for MULTIPLE Hopf, with user given coeff in aux.z; 
 z=aux.z; zl=length(z); ig=zeros(nu,1); 
 for j=1: zl; ig=ig+V(1:nu,j)*z(j)*exp(-2i*sign(om)*pi*t); end
 ig=ds*al*real(ig); 
end 

p.hopf.y=p.hopf.y+ig;hoplot(p,p.plot.pfig,p.plot.pcmp,p.hopf.aux); 
figure(p.plot.pfig); title('initial guess');
p.sw.para=para; p.sw.spcalc=0; p=hoMini(p); p.sol.ds=ds; par=p.u(p.nu+1:end); 
switch para
    case 3;  p.sol.restart=1;   % ---------- nat. param, use TOM 
  p.hopf.y=[p.hopf.y; 
            p.hopf.y(:,1)*ones(1,tl); % append aux vars
            p.hopf.T*ones(1,tl)]; % append T 
  fprintf('nat.parametr, lambda-guess=%g\n',lam+dlam*ds); 
  p.hopf.lam=lam+dlam*ds^2;       
  f=horhs(0,p.hopf.y(:,1),p,par); p.hopf.u0dot=f(1:nu); 
    case 4; p.sol.restart=0;           % --------- arcl, use honloopext 
  tau=zeros(1,nu*tl+2+nqh);     % the big tangent vector
  for i=1:tl; si=(i-1)*nu+1; % fill u-part of tangent 
      tau(si:si+nu-1)=ig(:,i)/ds;  end; % honorm(p,tau)
  tau(nu*tl+1+freeT)=dlam*ds^2; % the lam part, don't set T yet 
  tau=tau/honorm(p,tau); p.hopf.tau=tau; 
  lamg=lam+ds^2*tau(nu*tl+1+freeT); 
  fprintf('arc.parametr, lambda-guess=%g\n',lamg); 
  p.hopf.y=p.u(1:nu)*ones(1,tl)+p.sol.ds*v2tom(p,p.hopf.tau); 
  p=setlam(p,lamg); 
  p.hopf.y0d=sety0dot(p,p.hopf.y,par,p.hopf.T); % phase cond 
  if pcheck % check predictor 
    y0=p.u(1:nu)*ones(1,tl); % vals at HP (only u)
    y=y0+ds*v2tom(p,p.hopf.tau); %max(abs(p.hopf.tau))
    p.hopf.y=y; p.hopf.lam=lam+1*ds*tau(nu*tl+2); 
    hoplot(p,4,p.plot.pcmp,p.hopf.aux); title('initial guess'); 
    f=hoassempbc(p,p.hopf.y,p.hopf.T,p.hopf.lam); 
    normf=norm(f,inf); normy=max(max(abs(y(1:nu,:)-y0(1:nu,:)))); 
    fprintf('max|y-u0|=%g, res=%g\n',normy,normf); 
    p.hopf.y0d=sety0dot(p,p.hopf.y,par,p.hopf.T);
  end
  p.hopf.y=p.u(1:nu)*ones(1,tl); % reset y (only compute tangent for arclength!) 
   case 6;                     % --- Newton-Krylow shoting, after Sanchez-Net
  iT=aux.iT; p.hopf.nt=aux.nt; p.hopf.ntp=aux.ntp; p.nc.nq=aux.nq; 
  try; p.hopf.pskip=aux.pskip; catch p.hopf.pskip=10; end 
  try; p.hopf.sec=aux.sec; catch p.hopf.sec=0; end 
  try; p.hopf.ntmax=aux.ntmax; catch p.hopf.ntmax=aux.nt; end 
  try; p.hopf.al=aux.al; catch p.hopf.al=1; end 
  try; p.plot.view=aux.view; end 
  p.u(p.nu+iT)=p.hopf.T; p.hopf.iT=iT; par=p.u(p.nu+1:end); 
  p.hopf.G0=pderesi(p,[p.hopf.y(1:p.nu,1);par]); % for t-PCs 
  tau=[p.hopf.y(1:nu,1)-p.u(1:nu); zeros(2+p.nc.nq,1)]; 
  p.hopf.t=[];  p.hopf.y=[]; p.hopf.flcheck=0;  
  taun=xinormhoK(p,tau);  p.hopf.tau=tau/taun; 
  us=p.u; p.u(1:nu)=p.u(1:nu)+p.sol.ds*p.hopf.tau(1:nu); 
  mclf(1); hoplot(p,1,1,aux); 
  p.u=us; %plotsolu(p,p.hopf.tau,7,1,1); 
end 