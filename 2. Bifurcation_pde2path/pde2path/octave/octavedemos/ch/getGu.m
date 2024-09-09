function Gu=getGu(p,u,r)
% GETGU: compute jacobian of G, pde-part and auxiliart part separated.
if (p.sw.spcont==0) % regular continuation
   Gu=getGupde(p,u,r); % get pde part linearization
   if p.nc.nq>0 % derivatives w.r.t active auxiliary vars by finite diff.
     Gqw=zeros(p.nu+p.nc.nq,p.nc.nq); % (nu+nq) x nq
     for k=2:p.nc.nq+1 % derivative wrt parameters, skip primary parameter
       % ej(j,n)=jth canonical basis vector of length n
       r1=resi(p,u+p.nc.del*ej(p.nu+p.nc.ilam(k),length(u))); 
       Gqw(:,k-1)=(r1-r)/p.nc.del;
     end
   %  Gqw, pause 
     rq=r(p.nu+1:end); 
     if p.sw.qjac==1; qu=p.fuha.qfder(p,u); % analytical q_u 
     else; qu=zeros(p.nc.nq,p.nu); % numerical q_u 
         for j=1:p.nu 
            r1=p.fuha.qf(p,u+p.nc.del*ej(j,length(u))); 
            qu(:,j)=(r1-rq)/p.nc.del; 
         end
     end
     Gu=[[Gu; qu] Gqw];
   end 
   return;
end % regular continuation done.
if (p.sw.spcont==1); Gu=getGubpco(p,u,r);  end
if (p.sw.spcont==2); Gu=getGufoco(p,u,r);  end