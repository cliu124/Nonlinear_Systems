function p=dropp(p)
%remove phase invariance (modify p.mat.drop, p.mat.fill)
%drop the unknown corresponding to Im(u) at point p.pfn 
u=p.u(1:p.nu); par=p.u(p.nu+1:end); 
Drop_i=[1:p.nu-1]'; 
Drop_j=[1:p.nu/2+p.pfn-1 p.nu/2+p.pfn+1:p.nu]';
Drop_s=ones(p.nu-1,1); 
Drop=sparse(Drop_i,Drop_j,Drop_s,p.nu-1,p.nu); % nu-1 x nu matrix

p.mat.drop=Drop*p.mat.drop; % nu-1 x nu-1
Fill=transpose(Drop); p.mat.fill=p.mat.fill*Fill; p.nu=p.nu-1; 
p.u=zeros(p.nu+length(par),1);
p.u(1:p.nu)=Drop*u;  %currently not needed since always using dropp with zero u
p.u(p.nu+1:p.nu+length(par))=par; 
p=setfemops(p); 