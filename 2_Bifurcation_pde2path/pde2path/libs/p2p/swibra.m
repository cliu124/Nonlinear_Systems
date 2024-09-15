function q=swibra(dir,fname,varargin)
% swibra: branch switching (at simple BP) 
%
%  q=swibra(p)                         - bare syntax 
%  q=swibra(p,newdir,ds)               - set file.dir to newdir, use ds,      
%  q=swibra(p,newdir,ds,aux)           - add args in aux 
%  q=swibra(dir,fname)                 - load p from dir/fname
%  q=swibra(dir,fname,newdir)          - set file.dir to newdir
%  q=swibra(dir,fname,newdir,ds)       - also set initial stepsize to ds
%  q=swibra(dir,fname,newdir,ds,aux)   - use add.arg. from aux 
if ischar(dir); p=loadp(dir,fname);
else; p=dir;  % input is struct; shift other inputs; complicated due to backward comp. 
   if nargin>1; aux=varargin; 
      varargin{1}=fname; % newdir=fname in this case! 
      if ~isempty(aux); varargin{2}=aux{1}; end  
      if nargin>3; varargin{3}=aux{2}; end 
   end
end
try Xcont=p.sw.Xcont; catch Xcont=0; end 
if p.sol.ptype==3; fprintf('\n Hopf-point, use hoswibra\n'); q=[]; return; end; 
naux=length(varargin); aux=[]; 
if naux>1; p.sol.ds=varargin{2}; gotds=1; 
   if naux>2; aux=varargin{3}; end
else p.sol.ds=p.nc.dsmax/10; gotds=0; 
end
lam=getlam(p); 
try; fprintf(['lam=' num2str(lam), '; 4 smallest eigenvalues: ' num2str(p.sol.muv(1:4)) '\n']); 
catch; end
try mu2=aux.mu2; catch; mu2=p.nc.mu2; end; 
u=p.u; tau=p.tau; xi=p.sol.xi; lvar=length(varargin); 
aux=0; pstyle=p.plot.pstyle; hor=0; del=1e-3; % don't choose del too small!
if lvar>1; p.sol.ds=varargin{2}; gotds=1; lvar=lvar-1; 
   if lvar>1; aux=varargin{3}; end
else p.sol.ds=p.nc.dsmax/10; gotds=0; 
end 
if isfield(aux,'pstyle'); pstyle=aux.pstyle; end 
if isfield(aux,'hor'); hor=aux.hor; end 
u0d=tau(1:p.nu+p.nc.nq); al0=tau(p.nu+p.nc.nq+1); 
if hor==1;  fprintf('user requested horizontal swibra\n');
  tau1=[zeros(p.nu+p.nc.nq,1); 1]; 
else % calc phi1, al1, phi0 and then al0bar=a1 and al1bar
  r=resi(p,u); [Gu,Glam]=getder(p,u,r); % residual and jacobians 
  j=1; % use j=1 eigenval nearest to zero 
  if p.sw.eigsstart==1; vs=size(Gu,1); p.sw.evopts.v0=ones(vs,1)/vs; 
  else try; p.sw.evopts=rmfield(p.sw.evopts,'v0'); catch; end 
  end 
  M=getM(p);  % trivially extend M in case of auxiliary equations
  if (p.nc.nq>0); 
      if isfield(p.fuha,'qMfu'); [qL,qU,qD]=p.fuha.qMfu(p); else [qL,qU,qD]=stanqM(p); end 
      M=[[M(1:p.nu,1:p.nu) qU]; [qL qD]]; 
     % M=[[M zeros(p.nu,p.nc.nq)]; [zeros(p.nc.nq,p.nu) 0*speye(p.nc.nq)]];
  end   
  [phi1v,mu]=myeigs(Gu,M,j,0,p.sw.evopts,p); phi1=real(phi1v(:,j));
  fprintf('zero eigenvalue is %g\n',mu(j,j)); 
  fprintf('d/ds(lam)=%g\n',al0);
  [psi1v,mu]=myeigs(Gu',M,j,0,p.sw.evopts,p); 
  psi1=psi1v(:,j);psi1=psi1/(phi1'*psi1); 
  al1=psi1'*u0d;  phi0=(u0d-al1*phi1)/al0; 
  rp=resi(p,au2u(p,u2au(p,u)+del*phi1)); % residual and jacobian at u+del*phi1 
  [Gup,Glamp]=getder(p,au2u(p,u2au(p,u)+del*phi1),rp); Gud=Gup-Gu; 
  if any(Gud) 
    Glamd=Glamp-Glam; a1=psi1'*(Gud*phi1)/del; 
    b1=psi1'*(Gud*phi0+Glamd)/del; 
    al1b=-(a1*al1/al0+2*b1); 
    fprintf('al1=%g, a1=%g, b1=%g, al1b=%g\n',al1,a1,b1,al1b); 
    if al1b==0 fprintf('No distinct branch to switch to.');return; end
    tau1=[al1b*phi1+a1*phi0; a1];  
  else % trivial branch 
    tau1=[phi1; 0]; %fprintf('trivial swibra\n');
  end
end
tau1=tau1/xinorm(tau1,xi,p.nc.nq,p.sol.xiq); tau1=real(tau1); %tau1(end), pause
if Xcont>0; p.up=[tau1(1:p.nu);p.u(p.nu+1:end)]; end 
plotsolu(p,[tau1(1:p.nu);p.u(p.nu+1:end)],6,p.plot.pcmp,pstyle);
try; title(['\tau_' mat2str(p.plot.pcmp) ' at ' fname]); catch; end; 
if(p.sw.inter>1); p.sol.xi=asknu('xi',xi); 
    if ~gotds; p.sol.ds=asknu('ds',p.sol.ds); end 
end
p.tau=tau1; bs=p.branch(:,size(p.branch,2)); p=resetc(p); % set counters (note branch set below)
p.sol.ptype=-2; bs(1)=0; bs(2)=p.sol.ptype; % put last of p on new branch
if ~isempty(varargin); [p,ok]=setfn(p,varargin{1}); if ok~=1; q=p; return; end
else fprintf('warning: problem directory unchanged.\n'); end
p.branch=bs; p.sol.deta=0; p.file.count=0; 
p.sw.abs=1; % mark that we are right after branch-switching 
p.fuha.savefu(p); p.file.count=1;     
q=p;