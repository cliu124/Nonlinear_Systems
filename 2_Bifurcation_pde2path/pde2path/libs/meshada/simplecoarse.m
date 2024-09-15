function p=simplecoarse(p,varargin)
% simplecoarse: simple coarsening, 1D (hence OOPDE) 
%
% remove every 2nd mesh-point, unless the mesh-width is already>dxmax
aux=[]; if nargin>1; aux=varargin{1}; end 
onp=p.np; oldp=getpte(p);
u=p.mat.fill*p.u(1:p.nu); par=p.u(p.nu+1:end); 
tauf=p.mat.fill*p.tau(1:p.nu); taux=p.tau(p.nu+1:p.nu+p.nc.nq); % old full tau t
try dxmax=aux.dsmax; catch; try; dxmax=p.nc.dxmax; catch; dxmax=1; end; end 
dxv=diff(oldp); 
if mod(onp,2)==0;  rmc=2:2:onp; else; rmc=2:2:onp-1; end % removal candidates 
rmll=length(rmc); rml=[]; 
for i=1:rmll-1
  if dxv(rmc(i))<dxmax; rml=[rml rmc(i)]; end
end
keepl=setdiff(1:onp,rml); pon=oldp(keepl); np=length(keepl); 
trn=[1:np-1; 2:np]; % new 'triangulation' (trivial) 
p.pdeo.grid.p=pon; p.pdeo.grid.t=trn; p.np=np; p.pdeo.grid.e(1,2)=np; 
uc=zeros(p.nc.neq*np,1); 
for i=1:p.nc.neq % interpol to coarse mesh  
   uc((i-1)*np+1:i*np)=interp1(oldp(:),u((i-1)*onp+1:i*onp),pon);  
   taui((i-1)*np+1:i*np)=interp1(oldp(:),tauf((i-1)*onp+1:i*onp),pon); 
end
[p.mat.fill, p.mat.drop, p.nu]=getPerOp(p); 
p=oosetfemops(p); ucd=p.mat.drop*uc; p.u=[ucd; par]; p.nu=size(ucd,1); % coarse reduced u
p.tau=[p.mat.drop*taui';taux];  
p.r=resi(p,p.u); fprintf('inires=%g\n',norm(p.r,Inf));   
[p.u,r]=nloop(p,p.u); fprintf('res=%g\n',norm(r,p.sw.norm));
fprintf('coarsened from %i to %i\n',onp,p.np); 