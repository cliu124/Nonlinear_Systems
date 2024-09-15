% analytical 2nd Floquet multiplier for cGL with pBC (or NBC) 
function mvu=anafloq(rvek,s)
p=loadp('hom1d','hpt1'); par=getaux(p); muv=[]; 
for r=rvek; m=mfu(r,par,s); muv=[muv m]; end
end

function m=mfu(r,par,s)
nu=par(2); mu=par(3); c3=par(4); c5=par(5); 
a2=-c3/(2*c5)+s*sqrt(c3^2/(2*c5)^2+r);
om=nu-mu*a2; T=2*pi/om; hfak=r+3*a2-5*a2^2; m=exp(hfak*T); 
end