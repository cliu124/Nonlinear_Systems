function p=shinit(p,domsw,par,rad,nref,varargin) % Swift-Hohenberg as 2 component system
p=stanparam(p); p.nc.neq=2; p.nc.neig=4; % # evals to compute (few, for speed) 
p.fuha.outfu=@shbra; p.sol.ds=0.01; p.sol.dsmax=0.1; 
p.sw.bifcheck=2; p.sw.foldcheck=1; p.sw.spcalc=1; p.sw.sfem=-1;  
pde=stanpdeo2D(1,1,10); % dummy pdeo, next filled by 6-node triangulations 
switch domsw
  case 1; % disk 
  [p.nt,p.np,p.po,tri,p.efl,p.gfl]=trgl6_disk(nref,rad); % 6-node-triangulation 
  c3=c3fu(tri,p.nt); size(c3), c3=[c3, ones(size(c3,1),1)]; % convert to 3-node triangle-mesh 
  pde.grid.p=p.po'; pde.grid.t=c3'; p.pdeo=pde; % store 3-node-mesh (for plotting and Krot) 
  p.hofem.tri=tri; 
case 2;  % half disk    
 [p.nt,p.np,p.po,tri,p.efl,p.gfl]=trgl6_hdisk(nref,rad); % 6-node-triangulation 
 c3=c3fu(tri,p.nt); size(c3), c3=[c3, ones(size(c3,1),1)]; % convert to 3-node triangle-mesh 
 pde.grid.p=p.po'; pde.grid.t=c3'; p.pdeo=pde; % store 3-node-mesh (for plotting and Krot) 
 p.hofem.tri=tri; 
case 3;  % quarter disk    
 [p.nt,p.np,p.po,tri,p.efl,p.gfl]=trgl6_qdisk(nref,rad); % 6-node-triangulation 
 c3=c3fu(tri,p.nt); size(c3), c3=[c3, ones(size(c3,1),1)]; % convert to 3-node triangle-mesh 
 pde.grid.p=p.po'; pde.grid.t=c3'; p.pdeo=pde; % store 3-node-mesh (for plotting and Krot) 
 p.hofem.tri=tri; 
case 4; % disk sector
 aux=varargin{2}; p.pdeo=secpdeo2(rad,aux.nr,aux.nphi,aux.al,aux.phi,aux.jf);  
 p.np=p.pdeo.grid.nPoints; p.u=zeros(2*p.np,1); % prepare conversion to 6-nodes-triang. 
 p.t2sinterpol=0; % u re-initialized below, hence no need to interpolate u (or tau) to 6-nodes mesh, 
 p=tri2six(p.pdeo.grid.t',p.pdeo.grid.p',p); % convert to 6-node mesh 
end
p.plot.pstyle=2; p.plot.cm='parula'; p.plot.axis='image';  p.plot.bpcmp=11;
p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu; p.nc.lammin=-4; p.nc.lammax=2; 
u=0*ones(p.np,1); v=u; u0=[u v]; p.u=u0(:); 
p.u=[p.u; par']; p.nc.ilam=1; p.file.smod=5; p=setfemops(p); huclean(p); 
p.plot.auxdict={'\epsilon','\nu','q','sx','sy','srot'}; 
p.Om=sum(sum(p.mat.M(1:p.np,1:p.np),1)); % domain size, for normalization of L2-norm 