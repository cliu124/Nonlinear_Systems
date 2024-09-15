function [p,flag]=oomeshada(p,varargin)
% oomeshada: OO (and trullekul) mesh adaptation
%  p=oomeshada(p,varargin)
% 
% 1D,2D: classical adaption based on e2rs (element to refine selector), 
%        varargin: 'ngen',ngen  and 'sig',sig
%        or trullekul if p.sw.trul~=0
% 3D   : uniform refinement by OOPDE if p.sw.trul=0
%        otherwise trullekul, based on "pnorm metric" for z=p.trop.zfu(p) 
ngen=p.nc.ngen; noa=nargin-2; k=1;  imaxs=p.nc.imax; p.nc.imax=2*p.nc.imax; 
try trop=p.trop; catch trop=troptions3D();  end 
try sw=trop.sw; catch sw=0; end   % switch for tradapt 
try trsw=p.sw.trul; catch trsw=0;  end 
try maxt=p.nc.maxt; catch maxt=2*size(p.pdeo.grid.t,2); end 
try gosafe=p.sw.gosafe; catch; gosafe=0; end % safety-switch: if gosafe=1, successful 
% refinements are saved to disk to later reload the last successful one
% (necessary cause p.pdeo is a handle object and hence cannot be easily copied)
try swapfn=p.file.swapfn; catch; swapfn='maswap'; end 
while k<=noa 
   switch lower(varargin{k})       
       case 'ngen'; ngen=varargin{k+1}; k=k+2;
       case 'sig';  p.nc.sig=varargin{k+1}; k=k+2;
       otherwise; break
   end
end
oonp=p.np; flag=0; 
if gosafe; save(swapfn,'p'); end % save to later recover successful steps 
for j=1:ngen  % loop over refinements 
  u=p.u; uf=p.mat.fill*u(1:p.nu); oldp=p.pdeo.grid.p; onp=p.np; 
  par=u(p.nu+1:end); u=u(1:p.nu); tau=p.tau(1:p.nu); tauf=p.mat.fill*tau; 
  lamd=p.tau(p.nu+1:end);   
  switch size(oldp,1)
case 3; % 3D
  if trsw==0; p.pdeo.grid.refineMesh; % uniform refinement by OOPDE,  
  else % trullekrul 3D 
    [po,tr,ed]=getpte(p); xy=po'; tri=tr'; tri=tri(:,1:4); gr=p.pdeo.grid;    
    z=trop.zfu(p); eta=trop.etafu(p,onp); Nmetric=metric_pnorm(tri,xy,z,eta,trop.ppar); 
    bndmesh=[]; geomfunc=[]; bndmesh.fac=ed(1:3,:)'; bndmesh.IDs=ed(5,:)'; 
    bndmesh=bndmesh_polyhedron(tri,xy,bndmesh,trop);  
    if sw~=0; 
     [tri,xy,Nmetric,bndmesh,triQ]=tradapt(tri,xy,[Nmetric z],bndmesh,geomfunc,trop);
    else [tri,xy,Nmetric,bndmesh,triQ]=adapt_mesh(tri,xy,[Nmetric z],bndmesh,geomfunc,trop);
    end 
    nt=size(tri,1); ne=size(bndmesh.IDs,1); 
    gr.p=xy'; gr.t=[tri'; ones(1,nt)]; gr.e=bndmesh.fac'; 
    gr.e=[gr.e; zeros(1, ne)]; gr.e=[gr.e; 0*bndmesh.IDs']; 
    gr=trop.setids(gr); p.pdeo.grid=gr; 
   end
otherwise;  % 2D or 1D 
   if trsw~=0; % trullekrul 2D 
    [po,tr,ed]=getpte(p); xy=po'; tri=tr'; tri=tri(:,1:3); gr=p.pdeo.grid;       
    z=trop.zfu(p); eta=trop.etafu(p,onp); Nmetric=metric_pnorm(tri,xy,z,eta,trop.ppar); 
    bndmesh=[]; geomfunc=[]; bndmesh.tri=ed(1:2,:)'; bndmesh.IDs=ed(5,:)'; 
    bndmesh=bndmesh_polygon(tri,xy,bndmesh,trop);   
    if sw~=0; [tri,xy,Nmetric,bndmesh,triQ]=tradapt(tri,xy,[Nmetric z],bndmesh,geomfunc,trop);
    else [tri,xy,Nmetric,bndmesh,triQ]=adapt_mesh(tri,xy,[Nmetric z],bndmesh,geomfunc,trop);
    end 
    nt=size(tri,1); ne=size(bndmesh.IDs,1); gr.p=xy'; gr.t=[tri'; ones(1,nt)]; 
    gr.e=bndmesh.edg'; gr.e=[gr.e; zeros(2, ne)]; 
    gr.e=[gr.e; bndmesh.IDs']; gr=trop.setids(gr); % update boundary IDs 
    %gr.e, pause 
    p.pdeo.grid=gr;      
   else % OOPDE-error estimation and refinement      
    [p,idx]=p.fuha.e2rs(p,p.u); % select elements to refine 
    fprintf('refining %i elements\n',length(idx)); p.pdeo.grid.refineMesh(idx); 
     %p.pdeo.grid.refineMesh; 
   end   
  end
  npo=p.pdeo.grid.p; np=p.pdeo.grid.nPoints; p.sol.xi=1/np; ui=zeros(p.nc.neq*np,1); taui=ui;
  for i=1:p.nc.neq % interpolate u and tau to new mesh 
    switch size(oldp,1)  
      case 1; ui((i-1)*np+1:i*np)=interp1(oldp(:),uf((i-1)*onp+1:i*onp),npo);        
         taui((i-1)*np+1:i*np)=interp1(oldp(:),tauf((i-1)*onp+1:i*onp),p.pdeo.grid.p(:)); 
      case 2; xo=oldp(1,:); yo=oldp(2,:); 
         xn=p.pdeo.grid.p(1,:); yn=p.pdeo.grid.p(2,:); 
         ui((i-1)*np+1:i*np)=p2interpol(xn,yn,uf((i-1)*onp+1:i*onp),xo,yo,p); 
         taui((i-1)*np+1:i*np)=p2interpol(xn,yn,tauf((i-1)*onp+1:i*onp),xo,yo,p); 
      case 3; xo=oldp(1,:); yo=oldp(2,:); zo=oldp(3,:); 
         xn=p.pdeo.grid.p(1,:); yn=p.pdeo.grid.p(2,:); zn=p.pdeo.grid.p(3,:); 
         ui((i-1)*np+1:i*np)=p3interpol(xn,yn,zn,uf((i-1)*onp+1:i*onp),xo,yo,zo,p);        
         taui((i-1)*np+1:i*np)=p3interpol(xn,yn,zn,tauf((i-1)*onp+1:i*onp),xo,yo,zo,p);       
    end
  end 
  p.np=np; p.nu=p.nc.neq*p.np; p.u=[ui;par]; %size(p.u), size(taui), size(lamd)
  if size(lamd,1)>size(lamd,2); lamd=lamd'; end 
  p=box2per(p); p.tau=[p.mat.drop*taui;lamd']; 
  if size(oldp,1)==1 && p.sw.verb>1 % display new point distribution (1D) 
    figure(6); clf; plot(npo,'r*'); axis tight; 
    figure(7); clf; plot(npo,p.u(1:p.np),'*'); 
  else figure(6); clf; plotsol(p,6,1,1); axis(p.plot.axis); 
  end
  r=resi(p,p.u); fprintf('inires=%g\n',norm(r,Inf));   % was: r=p.r
  [p.u,r]=nloop(p,p.u); cres=norm(r,p.sw.norm); fprintf('res=%g\n',cres);
  if cres<p.nc.tol; flag=1; if gosafe; save(swapfn,'p'); end; end % refinement successful 
  nt=size(p.pdeo.grid.t,2); 
  if nt>maxt; fprintf('nt=%i>ntmax=%i, stopping adaption\n',nt,maxt);   break; end 
end
if gosafe; s=load(swapfn,'p'); p=s.p; end % recover last successful step 
p.nc.imax=imaxs; fprintf('refined from %i to %i\n',oonp,p.np);
