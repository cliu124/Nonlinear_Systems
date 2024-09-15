function p=meshadac(p,varargin)
% MESHADAC: mesh-adaption with coarsening
% 
%  p=meshadac(p,aux) 
% 
% aux=option-value pairs as in meshada
% 
% legacy setup: adapt mesh after coarsening, i.e. first interpolate to (coarse) 
% background-mesh (from setbmesh), then adapt. 
% in OOPDE: 1D: delete every 2nd mesh-point, unless mesh already too coarse (at
% places), i.e below p.nc.dsmax
% If p.sw.trul=1, then use trullekrul (2D and 3D), see tradapt. 
%
% See also meshada, setbmesh, oomeshadac
if p.sw.sfem<0; p=oomeshadac(p,varargin{:}); return; end % OO version
fprintf('   - current nt=%g, base nt=%g, ', p.mesh.nt, size(p.mesh.bt,2)); noa=nargin-1; 
nps=p.np; onp=size(p.mesh.bp,2); % save current points to later interpolate tau
fill_old = p.mat.fill; nu_old = p.nu;
lamd=p.tau(p.nu+p.nc.nq+1); % remember lam-direction 
upde=p.mat.fill*p.u(1:p.nu); % TD
uaux=p.u(p.nu+1:end); taux=p.tau(p.nu+1:p.nu+p.nc.nq); 
xo = p.mesh.p(1,:); yo = p.mesh.p(2,:);
xn = p.mesh.bp(1,:); yn = p.mesh.bp(2,:);
for i=1:p.nc.neq % interpol to coarse mesh 
    z=upde((i-1)*nps+1:i*nps); un=p2interpol(xn,yn,z,xo,yo); 
    uc((i-1)*onp+1:i*onp)=un; 
end 
p.mesh.p=p.mesh.bp; p.mesh.e=p.mesh.be; p.mesh.t=p.mesh.bt; p.np=onp; 
[p.mat.fill, p.mat.drop, p.nu]=getPerOp(p); 
uc=p.mat.drop*uc'; p.nu=size(uc,1); 
p.mesh.nt=size(p.mesh.t,2); p.u=[uc; uaux]; 
if(noa>0); p=meshada(p,varargin{:}); % refine the coarse mesh 
else p=meshada(p);
end
figure(p.plot.ifig);pdemesh(p.mesh.p,p.mesh.e,p.mesh.t);axis tight; 
plotsol(p,p.plot.pfig,p.plot.pcmp,p.plot.pstyle); xi=1/p.np; 
if(p.sw.inter>1); fprintf('old xi=%g, new xi=%g\n',p.sol.xi,xi); 
    xi=asknu('xi:',xi); end 
p.sol.xi=xi;
tau_old_full = fill_old*p.tau(1:nu_old);
tau_full=zeros(p.np*p.nc.neq,1); xn = p.mesh.p(1,:); yn = p.mesh.p(2,:);
for i=1:p.nc.neq % interpolate tau
  z=tau_old_full((i-1)*nps+1:i*nps); taun=p2interpol(xn,yn,z,xo,yo); 
  tau_full((i-1)*p.np+1:i*p.np) = taun; 
end
tau = p.mat.drop*tau_full; % new TD 
if p.nc.nq>0; tau(p.nu+1:p.nu+p.nc.nq)=taux; end % append aux variables to tau 
tau(p.nu+p.nc.nq+1)=lamd;  % append dlam coordinate and normalize
p.tau=tau/xinorm(tau,p.sol.xi,p.nu,p.sol.xiq); %p.fuha.headfu(p); 