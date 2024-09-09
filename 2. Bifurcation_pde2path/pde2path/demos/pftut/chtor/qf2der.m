function qu=qf2der(p,u)  
q1u=p.mat.vMs/p.vol; 
if p.nc.nq==1; qu=q1u; 
else uold=p.u(1:p.nu); uox=p.mat.Dphi*uold; q2u=uox'; qu=[q1u;q2u]; 
end
