function p=h2fdisk0(dir,pt,p,newdir)  % extend half-disk soln to full disk by 0. 
% Input p must contain mesh etc for full disk 
q=loadp(dir,pt,newdir); par=q.u(q.nu+1:end)';  % half-disk-soln
np=p.np; hp=q.pdeo.grid.p; nph=q.np; % half points
for i=1:np;
    mz=0; 
    po=p.pdeo.grid.p(:,i); x=po(1); y=po(2); 
    if x<0; x=-x; mz=1; end 
    dx=hp-repmat([x;y],1,nph); dxv=vecnorm(dx); 
    [mi,in]=min(dxv); % find point in p closest to point in q
    if mz==1; p.u(i)=0; p.u(np+i)=0; 
    else p.u(i)=q.u(in); p.u(np+i)=q.u(nph+in);         
    end
end
p.u(p.nu+1:end)=par; p.file=q.file; p.nc.dsmax=q.nc.dsmax; % cp stuff from half disk  
p=setfn(p,newdir); 