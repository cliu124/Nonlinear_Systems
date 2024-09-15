function plottauf(dir,sfname,wnr,cnr,pstyle)
% PLOTTAUF: plot tangent from dir/sfname.mat 
%
%  plottauf(dir,sfname,wnr,cnr,pstyle)
%
% wnr=window number, cnr=component number, pstyle=plot style
%
% See also plotsol, pdeplot, pdemesh
ffname=[dir '/' sfname '.mat'];tname=[dir sfname]; 
s=load(ffname,'p'); p=s.p; fprintf('lam=%g\n',getlam(p)); 
taupde=p.mat.fill*p.tau(1:p.nu);
figure(wnr); n0=(cnr-1)*p.np+1; n1=cnr*p.np;
if pstyle==1 pdemesh(p.mesh.p,p.mesh.e,p.mesh.t,taupde(n0:n1)); end 
if pstyle==2 pdeplot(p.mesh.p,p.mesh.e,p.mesh.t,'xydata',taupde(n0:n1)); end 
axis tight; box on; 
tname=[dir '/' sfname]; 
if(p.nc.neq>1) title(['\fontsize{16}\tau_' mat2str(cnr) ' at ' tname]); 
else title(['\fontsize{16}\tau at ' fname]); end 
if(p.plot.labelsw) xlabel('x','FontSize',p.plot.fs); 
  ylabel('y','FontSize',p.plot.fs); end;
set(gca,'FontSize',p.plot.fs); 
