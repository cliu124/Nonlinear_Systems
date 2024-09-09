function [h,hv,hk]=hfu(p,u) % harvesting values (and derivatives) for lvoc 
par=u(p.nu+1:end); al1=par(9); al2=par(10); 
v1=u(1); v2=u(p.np+1); k=lvcon(p,u); k1=k(1); k2=k(2); 
h1=v1^al1*k1^(1-al1); h2=v2^al2*k2^(1-al2); 
% pass back hv=[[h1v1 h1v2]; [h2v1 h2v2]] as vector [h1v1 h1v2 h2v1 h2v2]; 
h1v1=al1*v1^(al1-1)*k1^(1-al1); h1v2=0; h2v1=0; h2v2=al2*v2^(al2-1)*k2^(1-al2); 
h1k1=(1-al1)*v1^al1*k1^(-al1); h1k2=0; h2k1=0; h2k2=(1-al2)*v2^al2*k2^(-al2); 
h=[h1 h2]; hv=[h1v1 h1v2 h2v1 h2v2]; hk=[h1k1 h1k2 h2k1 h2k2]; 