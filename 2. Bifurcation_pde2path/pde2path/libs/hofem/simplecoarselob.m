function p=simplecoarselob(p)
% simplecoarselob: simple coarsening, Lobatto-nodes 
%
% remove every p.nc.redf -th mesh-point, unless the mesh-width is already>dxmax
try; redf=p.nc.redf; catch; redf=2; end 
onp=p.np; oldp=getpte(p); try; dxmax=p.nc.dxmax; catch; dxmax=1; end 
u=p.u(1:p.nu); par=p.u(p.nu+1:end); u=p.mat.fill*u; 
tauf=p.mat.fill*p.tau(1:p.nu); lamd=p.tau(end:end); % old full tau t
keepl=[1:redf:p.np-1 p.np]; %keepl, p.np
dxv=diff(oldp); 
for i=1:p.np-1    
  if dxv(i)>dxmax && ~ismember(i,keepl); keepl=[keepl i]; %i, keepl
  end
end
keepl=sort(keepl); npo=oldp(keepl); np=length(keepl); 
trn=[1:np-1; 2:np]; % new 'triangulation' (trivial) 
p.pdeo.grid.p=npo; p.pdeo.grid.t=trn; p.np=np; p.pdeo.grid.e(1,2)=np; 
p.hofem.xe=npo; p=two2lob(p); p=box2per(p);
np=p.np; npo=p.pdeo.grid.p; 
uc=zeros(p.nc.neq*np,1); taui=uc'; 
for i=1:p.nc.neq % interpol to coarse mesh  
   uc((i-1)*np+1:i*np)=interp1(oldp(:),u((i-1)*onp+1:i*onp),npo);  
   taui((i-1)*np+1:i*np)=interp1(oldp(:),tauf((i-1)*onp+1:i*onp),npo); 
end
if any(isnan(uc)) isnan(uc); end 
uc=p.mat.fill'*uc; 
p.u=[uc; par]; p.nu=size(p.mat.M,1); 
p.tau=[taui lamd]; % reduced u and tau 
fprintf('coarsened from %i to %i\n',onp,p.np); 