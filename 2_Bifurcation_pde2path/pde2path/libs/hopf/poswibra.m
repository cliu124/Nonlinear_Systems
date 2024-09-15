function p=poswibra(varargin)
% poswibra: branch-switching from periodic orbits (pitch/trans or period-doubling), 
% due to multipliers going through +-1
%
%  p=poswibra(dir,bptnr,newdir,ds)
%  p=poswibra(dir,bptnr,newdir,ds,aux) with 
%    aux.sw=1: enforce crit.multiplier at 1, 
%    aux.sw=-1: enforce crit.multiplier at -1 (PD) 
%    aux.coeff: case of higher-multipl. crit. multiplier; coeff for crit vectors
aux=[]; 
if ischar(varargin{1}) % dir/pt syntax
   dir=varargin{1}; fname=varargin{2}; ndir=varargin{3}; ds=varargin{4}; 
   if nargin>4; aux=varargin{5}; end 
   p=loadp(dir,fname);    
else p=varargin{1}; ndir=varargin{2};  ds=varargin{3};  
   if nargin>3; aux=varargin{4}; end  
end 
lam=getlam(p); fprintf('poswibra, lambda=%g\n',lam); 
dlam=0; if isfield(aux,'dlam'); dlam=aux.dlam; end 
sw=0; if isfield(aux,'sw'); sw=aux.sw; end 
c=1; if isfield(aux,'coeff'); c=aux.coeff; end 
try; nqh2=length(p.hopf.ilam); catch; nqh2=0; end 
ho=p.hopf; yorg=ho.y; 
lam=getlam(p); [f,jac]=tomassempbc(p,p.hopf.y,p.hopf.T,lam); 
[muv1,muv2,ind,mo,Vs1, Vs2]=floq(p,jac); 
try mu1=muv1(end-1); catch mu1=-999; end 
try mu2=muv2(1); catch mu2=999; end 
[m1,idx]=min([abs(mu1-1),abs(mu2-1),abs(mu1+1),abs(mu2+1)]);
if sw==1;  [m1,idx]=min([abs(mu1-1),abs(mu2-1)]); % enforce mu_c=1 
elseif sw==-1; [m1,idx]=min([abs(mu1+1),abs(mu2+1)]); end  % enforce mu_c=-1 
%sw, idx, pause 
switch idx
    case {1,3}; mu=mu1; y0=Vs1(:,end-1);        
    case {2,4}; mu=mu2; y0=zeros(p.nu,1); for i=1:length(c); y0=y0+c(i)*Vs2(:,i); end 
end
y0=real(y0); 
if real(mu)>0 % propagate y0 via Jac! 
    ytau=pobifpred(p,mu,y0,jac); % norm(ytau,'inf'), ds
    ynew=yorg+ds*ytau; ho.y=ynew; p.hopf.y=ynew; 
    tl=ho.tl; nu=p.nu; nqh=ho.nqh; 
    tau=zeros(1,nu*tl+2+nqh);     % the big tangent vector
    for i=1:tl; si=(i-1)*nu+1; % fill u-part of tangent 
          tau(si:si+nu-1)=ds*ytau(:,i);  end; % honorm(p,tau)
    tau(nu*tl+2)=dlam; % the lam part, don't set T yet 
else % PD-case: double tv, T, yorg etc 
  ho.T=2*ho.T; ho.t=[ho.t(1:end)/2, 0.5+ho.t(2:end)/2]; ho.tl=2*ho.tl-1; 
  yorg=[yorg(:,1:end-1) yorg]; 
  ytau1=pobifpred(p,mu,y0,jac); % propagate y0 for half a (new) period 
  y1=-y0; ytau2=pobifpred(p,mu,y1,jac); % other half period
  ytau=[ytau1(:,1:end-1) ytau2]; 
  ynew=yorg+1*ds*ytau; ho.y=ynew; 
  tl=ho.tl; nu=p.nu; nqh=ho.nqh; 
  tau=zeros(1,nu*tl+2+nqh2);     % the big tangent vector
  for i=1:tl; si=(i-1)*nu+1; % fill u-part of tangent 
          tau(si:si+nu-1)=ds*ytau(:,i);  end; % honorm(p,tau)
  tau(nu*tl+2)=dlam; % the lam part, don't set T yet   
end
p.hopf=ho; hoplot(p,10,1); p=resetc(p); p=setfn(p,ndir); 
p.hopf.tau=tau; par=p.u(p.nu+1:end); 
lamg=p.hopf.lam+dlam; p=setlam(p,lamg); p.sol.ds=ds; 
p.hopf.y0d=sety0dot(p,p.hopf.y,par,p.hopf.T); % phase cond 
%[pc,pc_y]=hopc(p,p.hopf.y,p.hopf.T,lamg); pc, norm(pc_y,'inf'),  
p.hopf.y=yorg; p.branch=[bradat(p); p.fuha.outfu(p,p.u)];  
p.fuha.savefu(p); p.file.count=1; p.hopf.indini=0; 
%p.hopf.y=ynew; 
end