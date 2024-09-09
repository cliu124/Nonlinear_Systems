function Gu=getGupde(p,u,r)
% convenience function to get jacobian for pde part used in getGu
% overload here for "phase-fix"
global pj; % for call of numjac
if(p.sw.jac==1) % pde-part analytically 
  upde=p.mat.fill*p.u(1:p.nu); 
  if p.sw.sfem~=1 % use full jac
    bc=p.fuha.bcjac(p,u); 
    [cj,aj,bj]=p.fuha.Gjac(p,u); zerov=zeros(p.nc.neq,1); 
    [Gu,dum]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,cj,aj,zerov,upde); 
    if(any(bj)); Kadv=assemadv(p.mesh.p,p.mesh.t,bj); Gu=Gu-Kadv; end 
    Gu=filltrafo(p,Gu);
  else % use "simple jac". p.mat.fill is built in here already
      Gu=p.fuha.sGjac(p,u);
  end
else 
    bc=p.fuha.bcjac(p,u);
    %[S,dum]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,...
   %         0,ones(p.nc.neq*p.nc.neq,1),ones(p.nc.neq,1),u(1:p.nu));
    %S=filltrafo(p,S);
    M0=p.mat.M(1:p.nu/2,1:p.nu/2)>0; S=[[M0 M0]; [M0 M0]]; 
    thresh=p.nc.del*ones(p.nu,1); pj=p; pj.u=u;
    [Gu,dum,dum]=numjac('resinj',0,u(1:p.nu),r(1:p.nu),thresh,[],0,S,[]); 
end
if p.pffac>0; % phase-fix at point p.ppn
    Gu(p.nu/2+p.pfn,:)=zeros(1,p.nu);
    Gu(p.nu/2+p.pfn,p.nu/2+p.pfn)=p.pffac;
end


