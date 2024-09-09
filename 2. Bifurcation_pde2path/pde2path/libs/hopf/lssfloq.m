function x=lssfloq(AA,BB,b,p) 
% lssfloq: AMG interface for it. soln of LSS for Floq-comp (eperimental!) 
n=p.nu; m=size(AA,3)+1; x=b; 
for i=1:m-1; 
   x=AA(:,:,i)*x; 
   [x,p]=lssAMG(sparse(BB(:,:,i)),x,p); 
   %x=BB(:,:,i)\x; 
end
    