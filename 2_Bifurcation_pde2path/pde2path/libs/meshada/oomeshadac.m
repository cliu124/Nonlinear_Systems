function p=oomeshadac(p,varargin)
% OOMESHADAC: adapt mesh after coarsening
%  original: first interpolate to (coarse) standard-mesh (from setbmesh), then adapt. 
%  now also with simplecoarse (1D), or trullekrul (if p.sw.trul=1) 
%
%  p=oomeshadac(p,varargin) 
% meshadapt from background mesh
%
% See also oomeshada, setbmesh, 
try; trul=p.sw.trul; catch trul=0; end 
if trul; p=oomeshadactr(p,varargin); return; end; % use trullekrul
try; scoarse=p.sw.scoarse; catch; scoarse=1; end % simple coarsening 
noa=nargin-1;  [po,t,e]=getpte(p);  
onp=size(po,2); oldp=po; % save current points to later interpolate tau
lamd=p.tau(p.nu+p.nc.nq+1); % remember lam-direction 
u=p.mat.fill*p.u(1:p.nu); tau=p.mat.fill*p.tau(1:p.nu); % old full tau 
uaux=p.u(p.nu+1:end); taux=p.tau(p.nu+1:p.nu+p.nc.nq); 
if scoarse && size(po,1)==1 
  p=simplecoarse(p); 
else  
    fprintf('   - current nt=%g, base nt=%g, ', size(t,2), size(p.mesh.bt,2));         
    p.pdeo.grid.p=p.mesh.bp; p.pdeo.grid.t=p.mesh.bt; % set mesh to background
    p.pdeo.grid.e=p.mesh.be; p.np=p.pdeo.grid.nPoints; np=p.np; 
    uc=zeros(p.nc.neq*np,1); 
    for i=1:p.nc.neq % interpol to coarse mesh  
      switch size(oldp,1) %p.ndim
          case 1; uc((i-1)*np+1:i*np)=interp1(oldp(:),u((i-1)*onp+1:i*onp),p.mesh.bp);        
          case 2; xo=oldp(1,:); yo=oldp(2,:); xn=p.mesh.bp(1,:); yn=p.mesh.bp(2,:);
             uc((i-1)*np+1:i*np)=p2interpol(xn,yn,u((i-1)*onp+1:i*onp),xo,yo); 
      end
    end 
    [p.mat.fill, p.mat.drop, p.nu]=getPerOp(p); 
    p=oosetfemops(p); ucd=p.mat.drop*uc; p.u=[ucd; uaux]; p.nu=size(ucd,1); % coarse reduced u 
end 
if(noa>0); [p,flag]=oomeshada(p,varargin{:}); % refine the coarse mesh 
else [p,flag]=oomeshada(p); end
if flag>0 
np=p.np; taun=zeros(p.nc.neq*np,1); [npo,t,e]=getpte(p); 
for i=1:p.nc.neq % interpolate tau
  switch size(oldp,1) 
      case 1; newtau=interp1(oldp(:),tau((i-1)*onp+1:i*onp),npo(:));  
          taun((i-1)*np+1:i*np)=newtau; 
      case 2; xo=oldp(1,:); yo=oldp(2,:); 
         xn=npo(1,:); yn=npo(2,:);
         taun((i-1)*np+1:i*np)=p2interpol(xn,yn,tau((i-1)*onp+1:i*onp),xo,yo); 
  end
end
tau=p.mat.drop*taun; % reduced tau (for pBC) 
if p.nc.nq>0; tau(p.nu+1:p.nu+p.nc.nq)=taux; end % append aux variables to tau 
tau(p.nu+p.nc.nq+1)=lamd;  % append dlam coordinate and normalize
p.tau=tau/xinorm(tau,p.sol.xi,p.nu,p.sol.xiq);%p.fuha.headfu(p); 
end