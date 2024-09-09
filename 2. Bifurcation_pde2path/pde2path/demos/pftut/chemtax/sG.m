function r=sG(p,u) % pde for chemotaxis model 
% u_t=0.25*Lap u-lam*div(u*grad v)+r*u(1-u),   v_t=Lap v+(u/(1+u)-v) 
lam=u(p.nu+1:end);d=0.25;r=1.52;u=u(1:p.nu);n=p.np;u1=u(1:n); u2=u(n+1:2*n); 
f1=r*u1.*(1-u1); f2=u1./(1+u1)-u2; f=[f1;f2]; % semilin.nonlinearity 
ut=p.mat.p2c*u1; 
gr=p.pdeo.grid; fem=p.pdeo.fem; 
cc=p.fuha.cfu(ut,lam); % quasilin. coefficient, here simply c(u)=u
[K12,~,~]=fem.assema(gr,cc,0,0); % assemble matrix for cross-diff 
K=p.mat.K; r=[d*K -lam*K12; 0*K K]*u-p.mat.M*f; % putting rhs together