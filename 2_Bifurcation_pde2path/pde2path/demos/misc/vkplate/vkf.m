function [c,a,f,b]=vkf(p,u) 
% von karman-plate as a 10 component system including some regularization
% i.e., aux variables like u5 defined as u5=(1-del Lap)^(-1)\pa_x^2 u1 
lam=u(p.nu+1); del=u(p.nu+2); % regularization parameters of aux varia
c=zeros(400,1); 
% first column of C as written in manual 
c(1)=1; c(4)=1; c(5)=lam; c(17)=1; c(24)=1; c(26)=0.5;c(27)=0.5; 
c(45)=1; c(48)=1; % second column 
c(89)=1; c(92)=1; c(109)=1; c(116)=1; c(118)=0.5; c(119)=0.5; % third
c(133)=1; c(136)=1;     % 4th
c(177)=del; c(180)=del; % 5th 
c(221)=del; c(224)=del; 
c(265)=del; c(268)=del; 
c(309)=del; c(312)=del; 
c(353)=del; c(356)=del; 
c(397)=del; c(400)=del;   % to 10th 
a=zeros(100,1); a(11)=1; a(33)=1; a(45)=1; a(56)=1; 
a(67)=1; a(78)=1; a(89)=1; a(100)=1;
b=0; 
u5=pdeintrp(p.mesh.p,p.mesh.t,u(4*p.np+1:5*p.np)); 
u6=pdeintrp(p.mesh.p,p.mesh.t,u(5*p.np+1:6*p.np)); 
u7=pdeintrp(p.mesh.p,p.mesh.t,u(6*p.np+1:7*p.np)); 
u8=pdeintrp(p.mesh.p,p.mesh.t,u(7*p.np+1:8*p.np)); 
u9=pdeintrp(p.mesh.p,p.mesh.t,u(8*p.np+1:9*p.np)); 
u10=pdeintrp(p.mesh.p,p.mesh.t,u(9*p.np+1:10*p.np)); 
f1=-(u5.*u9-2*u7.*u10+u6.*u8); 
f2=u5.*u6-u7.^2;  zv=zeros(1,p.mesh.nt); 
f=[zv;f1;zv;f2;zv;zv;zv;zv;zv;zv];