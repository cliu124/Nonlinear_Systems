function tintplot1d(dir,si,incr,nt,wnr,cmp,vv)
% tintplot1d: plot output from time integration 
% 
%  tintplot1d(dir,si,incr,nt,wnr,cmp,vv)
%
% dir=out-dir of tint, si=startindex, incr=increment
% nt=# t levels, wnr=window-nr, cmp=sol-component, vv=view
p=loadp(dir,'pt0'); nu=p.nu; sol.x=zeros(1,nt); 
sol.y=zeros(nu,nt); 
for i=1:nt;    
    pt=[dir '/pt' mat2str(si+(i-1)*incr)]; 
    tmp=load(pt); tmpp=tmp.p; %plotsol(p); pt, pause
    sol.x(i)=tmpp.t; sol.y(:,i)=tmpp.u(1:p.nu); 
end
xtplot(p,sol,wnr,cmp,vv,''); 
end 