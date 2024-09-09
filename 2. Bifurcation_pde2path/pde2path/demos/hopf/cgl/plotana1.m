% plot analytical Hopf-branch for cGL for comparison with numerics 
function plotana1(k2,rstep,ps,sw,ms)
p=loadp('hom1d','hpt1'); par=getaux(p); 
rv1=k2:-rstep:k2-0.25; av1=0*rv1; Tv1=av1; rl=length(rv1); 
for i=1:rl; r=rv1(i); [a,T]=afu(r,par,k2,-1); av1(i)=a; Tv1(i)=T; end
rv2=k2-0.25:rstep:k2+2; av2=0*rv2; Tv2=av2; rl=length(rv2); 
for i=1:rl; r=rv2(i); [a,T]=afu(r,par,k2,1); av2(i)=a; Tv2(i)=T; end
rv=[rv1, rv2]; Tv=[Tv1 Tv2]; av=[av1 av2]; 
set(0,'DefaultLineLineWidth', 2); figure(3); 
if sw==1; h=plot(rv,av,ps); 
    else h=plot(rv,Tv,ps); end 
    set(h, 'Markersize',ms);
set(0,'DefaultLineLineWidth', 1);
end

function [a,T]=afu(r,par,k2,s)
nu=par(2); mu=par(3); c3=par(4); c5=par(5); 
a2=-c3/(2*c5)+s*sqrt(c3^2/(2*c5)^2+r-k2);
a=sqrt(a2); om=nu-mu*a2; T=2*pi/om; 
end