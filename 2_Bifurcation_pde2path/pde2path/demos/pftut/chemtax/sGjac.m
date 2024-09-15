function Gu=sGjac(p,u) % chemotaxis model 
% u_t=0.25*Lap u-lam*div(u*grad v)+r*u(1-u), v_t=Lap v+(u/(1+u)-v) 
d=0.25; r=1.52; lam=p.u(p.nu+1); u=u(1:p.nu); n=p.np; u1=u(1:n); u2=u(n+1:2*n); 
f1u=r*(1-2*u1); f1v=0*f1u; f2u=1./((1+u1).^2); f2v=-ones(p.np,1); 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)]; % Jac of semilin. nonlin. 
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
ut=p.mat.p2c*u1; gr=p.pdeo.grid; fem=p.pdeo.fem; cc=p.fuha.cfu(ut,lam); 
K=p.mat.K; [K12,~,~]=fem.assema(gr,cc,0,0);  % cross-diffusion 
jacsw=0; 
switch jacsw; 
    case 0 % approximate way, does not work for all branches here 
  vx=p.mat.Dx*u2; vy=p.mat.Dy*u2; % 1st derivatives as coefficients   
  KK=p.mat.Kx*spdiags(vx,0,n,n)+p.mat.Ky*spdiags(vy,0,n,n); 
  Gu=[d*K+lam*KK, -lam*K12;  0*K, K]-p.mat.M*Fu; 
    case 1 % using getKuvd to obtain \pa_u \div(c(u)\nabla v) via numjac
  Kuvd=getKuvd(p,lam,u1,u2); 
  Gu=[d*K-lam*Kuvd -lam*K12;  0*K K]-p.mat.M*Fu;  
end