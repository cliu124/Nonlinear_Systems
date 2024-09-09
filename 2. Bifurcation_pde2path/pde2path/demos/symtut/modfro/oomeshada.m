%%
% oomeshada: OO mesh adaption
%  p=oomeshada(p,varargin)
%
% here adapted to take care of the u0 adaption
function p=oomeshada(p,varargin)
ngen=p.nc.ngen; noa=nargin-1; k=1;  
while k<=noa 
   switch lower(varargin{k})       
       case 'ngen'; ngen=varargin{k+1}; k=k+2;
       case 'sig';  p.nc.sig=varargin{k+1}; k=k+2;
       otherwise; break
   end
end
for j=1:ngen  % loop over refinements 
  u=p.u; uf=p.mat.fill*u(1:p.nu); oldp=p.pdeo.grid.p; onp=p.np; 
  par=u(p.nu+1:end); u=u(1:p.nu); tau=p.tau(1:p.nu); tauf=p.mat.fill*tau; 
  lamd=p.tau(p.nu+1:end);   
  [p,idx]=p.fuha.e2rs(p,p.u); % select elements to refine 
  fprintf('refining %i elements\n',length(idx)); 
  p.pdeo.grid.refineMesh(idx); npo=p.pdeo.grid.p; np=p.pdeo.grid.nPoints; 
  ui=zeros(p.nc.neq*np,1); taui=ui; 
  for i=1:p.nc.neq % interp. u and tau to new mesh (full mesh) 
    switch size(oldp,1) %p.ndim
      case 1; ui((i-1)*np+1:i*np)=interp1(oldp(:),uf((i-1)*onp+1:i*onp),npo);        
         taui((i-1)*np+1:i*np)=interp1(oldp(:),tauf((i-1)*onp+1:i*onp),p.pdeo.grid.p(:)); 
      case 2; xo=oldp(1,:); yo=oldp(2,:); 
         xn=p.pdeo.grid.p(1,:); yn=p.pdeo.grid.p(2,:); 
         ui((i-1)*np+1:i*np)=p2interpol(xn,yn,uf((i-1)*onp+1:i*onp),xo,yo); 
         taui((i-1)*np+1:i*np)=p2interpol(xn,yn,tauf((i-1)*onp+1:i*onp),xo,yo); 
    end
  end 
  p.np=np; p.nu=p.nc.neq*p.np; p.u=[ui;par]; 
  p=box2per(p); p.tau=[p.mat.drop*taui;lamd]; 
  if size(oldp,1)==1 && p.sw.verb>1 % check new point distribution
    figure(6); clf; plot(npo,'r*'); axis tight; 
    figure(7); clf; plot(npo,p.u(1:p.np),'*'); 
  else figure(6); clf; plotsol(p,6,1,1); axis(p.plot.axis); 
  end
  p.r=resi(p,p.u); fprintf('inires=%g\n',norm(p.r,Inf)); 
  [p.u,r]=nloop(p,p.u); fprintf('res=%g\n',norm(r,p.sw.norm));
end