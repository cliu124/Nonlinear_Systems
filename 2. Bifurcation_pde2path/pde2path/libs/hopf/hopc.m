function [pc,pc_y]=hopc(p,y,T,lam) 
% hopc: phase condition
nu=p.nu+p.nc.nq; withM=0; 
h=p.hopf.t(2:end)-p.hopf.t(1:end-1); os=1; % non-uniform grid
if withM; pci=sum((p.mat.M*p.hopf.y0d(1:nu,1:p.hopf.tl-os)).*y(1:nu,1:p.hopf.tl-os),1); 
else  pci=sum(p.hopf.y0d(1:nu,1:p.hopf.tl-os).*y(1:nu,1:p.hopf.tl-os),1);
    % sum along first components
end 
pc=sum(h(1:end+1-os).*pci);  
pc_y=zeros(1,nu*p.hopf.tl); 
for i=1:p.hopf.tl-os; 
    ic=(i-1)*nu; 
    if withM; pc_y(ic+1:ic+nu)=h(i)*(p.mat.M*p.hopf.y0d(1:nu,i))'; 
    else pc_y(ic+1:ic+nu)=h(i)*(p.hopf.y0d(1:nu,i))'; 
    end
end
pcfac=10; if isfield(p.hopf,'pcfac'); pcfac=p.hopf.pcfac; end % compatibility
pc=pcfac*pc; pc_y=pcfac*pc_y; 
return
% some testing .. 
tv=p.hopf.t; tl=length(tv); zv=zeros(2,tl); 
for i=1:tl; 
    zv(1,i)=norm(p.hopf.y0d(1:nu,i)); 
    zv(2,i)=(p.hopf.y0d(1:nu,i)'*y(1:nu,i)); 
end 
figure(101); clf; plot(tv,zv(1,:), tv, zv(2,:)); legend('|u_t|','u_t*u');
figure(102); clf; plot(tv(1:end-1),pci); legend('pci'); 
pc, norm(y(1:nu,1)-y(1:nu,end))
drawnow, pause
end 

