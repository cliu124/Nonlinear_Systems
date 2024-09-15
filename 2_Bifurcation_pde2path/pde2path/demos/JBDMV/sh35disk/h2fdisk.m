function p=h2fdisk(dir,pt,p,newdir)  % mirror half-disk to full 
% input p must contain mesh etc for full disk 
q=loadp(dir,pt,newdir); par=q.u(q.nu+1:end)';  % half-disk-soln
np=p.np; hp=q.pdeo.grid.p; nph=q.np; % half points
for i=1:np;
    po=p.pdeo.grid.p(:,i); x=po(1); y=po(2); 
    if x<0; x=-x; end 
    dx=hp-repmat([x;y],1,nph); dxv=vecnorm(dx); 
    [mi,in]=min(dxv); % find point in p closest to point in q
    p.u(i)=q.u(in); p.u(np+i)=q.u(nph+in); 
end
p.u(p.nu+1:end)=par; p.file=q.file; p.nc.dsmax=q.nc.dsmax; % cp stuff from half disk  
p=setfn(p,newdir); 