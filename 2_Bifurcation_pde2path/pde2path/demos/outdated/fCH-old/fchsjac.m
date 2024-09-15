function Gu=fchsjac(p,u)  
par=u(p.nu+1:end); u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
eta1=par(1); ga=par(2); eps=par(3); eta2=par(5); etad=eta1-eta2; 
[w,wp,wpp,wppp]=p.fuha.wfu(u1,p); ov=ones(p.np,1); 
f1u=-wpp; f1v=-ov; 
f2u=-wppp.*u2+eps*etad*wpp; 
f2v=-wpp+eps*eta1*ov; 
Fu=[[spdiags(f1u,0,p.np,p.np),spdiags(f1v,0,p.np,p.np)];
    [spdiags(f2u,0,p.np,p.np),spdiags(f2v,0,p.np,p.np)]]; 
Gu=eps^2*p.mat.K-p.mat.M*Fu; 
end
