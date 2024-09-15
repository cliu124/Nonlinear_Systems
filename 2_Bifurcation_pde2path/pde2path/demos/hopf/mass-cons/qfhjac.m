function qjac=qfhjac(p,y) % u-derivatives of qfh
qjac=zeros(1,p.hopf.tl*p.nu); j=p.mat.M(1:p.np,1:p.np)*ones(p.np,1)/p.vol; 
for i=1:p.hopf.tl;  % same derivative at each time slice 
    qjac((i-1)*p.nu+1:i*p.nu)=([j; j])'; 
end 