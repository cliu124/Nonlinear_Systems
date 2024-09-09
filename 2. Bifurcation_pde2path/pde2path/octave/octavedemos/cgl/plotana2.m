% plot analytic sol, here general parameter dependence par from mi to ma
function plotana2(p,k2,ilam,mi,ma,st,ps,sw,s) 
par=getaux(p); pv=mi:st:ma; pvl=length(pv); av=0*pv; Tv=av; 
for i=1:pvl; par(ilam)=pv(i); 
    [a,T]=afu(par,k2,s); av(i)=a; Tv(i)=T; end
set(0,'DefaultLineLineWidth', 2); figure(3); 
if sw==1; plot(pv,av,ps); else plot(pv,Tv,ps); end 
set(0,'DefaultLineLineWidth', 1);
end

function [a,T]=afu(par,k2,s)
r=par(1); nu=par(2); mu=par(3); c3=par(4); c5=par(5); 
a2=-c3/(2*c5)+s*sqrt(c3^2/(2*c5)^2+r-k2); 
a=sqrt(a2); om=nu-mu*a2; T=2*pi/om; 
end