function [cj,aj,bj]=vkjac(p,u) 
% jacobian for von karman-plate
[cj,a,f,b]=vkf(p,u);  % re-use c from vkf
bj=0;fu=zeros(p.nc.neq^2,p.mesh.nt); v1=ones(1,p.mesh.nt); 
fu(11,:)=-v1; fu(33,:)=-v1; fu(45,:)=-v1; fu(56,:)=-v1; 
fu(67,:)=-v1; fu(78,:)=-v1; fu(89,:)=-v1; fu(100,:)=-v1;
u5=pdeintrp(p.mesh.p,p.mesh.t,u(4*p.np+1:5*p.np)); 
u6=pdeintrp(p.mesh.p,p.mesh.t,u(5*p.np+1:6*p.np)); 
u7=pdeintrp(p.mesh.p,p.mesh.t,u(6*p.np+1:7*p.np)); 
u8=pdeintrp(p.mesh.p,p.mesh.t,u(7*p.np+1:8*p.np)); 
u9=pdeintrp(p.mesh.p,p.mesh.t,u(8*p.np+1:9*p.np)); 
u10=pdeintrp(p.mesh.p,p.mesh.t,u(9*p.np+1:10*p.np));
fu(42,:)=-u9; fu(52,:)=-u8; fu(62,:)=2*u10; 
fu(72,:)=-u6; fu(82,:)=-u5; fu(92,:)=2*u7;
fu(44,:)=u6; fu(54,:)=u5; fu(64,:)=-2*u7; 
aj=-fu; 