function p=twswibra(dir,fname,spar,kwnr,varargin)
% TWSWIBRA: branch switching from Hopf points to TWs in O(2) symmetry
% 
%  p=twswibra(dir,fname,speedpar,kwnr,newdir,aux)
%  
%  speedpar=index of speed s
%  kwnr=spatial wave-number (s.t. period T=p.L/(kwnr*p.u(p.nu+p.spar))
%  varargin=aux may include: 
%    z for multiple Hopf: coefficients to mix EVecs at bif 

if nargin>5; aux=varargin{2}; ndir=varargin{1}; else; aux=0; end 
p=loadp(dir,fname); nu=p.nu;  
if(~isempty(varargin)); [p,ok]=setfn(p,varargin{1}); 
    if ~ok; fprintf('could not change dir-name!\n'); return; end
else fprintf('warning: problem directory unchanged.\n'); end
lam=getlam(p); p=resetc(p); brout=[bradat(p); p.fuha.outfu(p,p.u)]; brout(1,1)=0; 
p.branch=brout; p.file.count=0; p.sol.ptype=2; p.fuha.savefu(p); p.file.count=1; 
r=resi(p,p.u); Gu=getGu(p,p.u,r); [ineg,muv,V]=spcalc(Gu,p,p.sol.j0); 
evnr=1; om=imag(muv(evnr)); mu0=real(muv(evnr)); 
idx1=abs(real(muv))<p.nc.mu2; muv=muv(idx1); V=V(:,idx1); % extract central evals 
[mus,idx]=sort(abs(imag(muv)-om));  % sort those close to om_0 to front 
muv=muv(idx); V=V(:,idx); phi=V(:,evnr); phi=phi/norm(phi,2); 
s=abs(om)/kwnr; % speed 
fprintf('real(mu)=%g, imag(mu)=%g\n',mu0,imag(muv(1))); fprintf('lam=%g, om=%g, speed=%g \n',lam, om,s); 
if ~isfield(aux,'z') ig=real(phi(1:nu));  % standard initial guess
else % initial guess, version for MULTIPLE Hopf, with user given coeff in aux.z; 
 z=aux.z; zl=length(z); ig=zeros(nu,1); 
 for j=1: zl; ig=ig+V(1:nu,j)*z(j); end
 ig=real(ig); 
end 
tau1=[ig; 0;0]/xinorm(ig,p.sol.xi,p.nc.nq,p.sol.xiq); % tangent 
plotsolu(p,[tau1;p.u(p.nu+1:end)],6,p.plot.pcmp,p.plot.pstyle); 
try; title(['\tau_' mat2str(p.plot.pcmp) ' at ' fname]); catch; end; 
p.tau=tau1; bs=p.branch(:,size(p.branch,2)); p=resetc(p); p.sol.restart=0; % set counters 
p.sol.ptype=-2; bs(1)=0; bs(2)=p.sol.ptype; % put last of p on new branch
p.branch=bs; p.sol.deta=0; p.u(p.nu+spar)=s; p.kwnr=kwnr; p.spar=spar; 
p.fuha.outfu=@hobratw; % computes T from k,s and domain size (1D) 
p.fuha.savefu(p); p.file.count=1; 
