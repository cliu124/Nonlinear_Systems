function Gu=sGjac(p,u) % Jac for MWBS'97 semiconductor model 
n=p.np; [f1u,f1v,f2u,f2v]=nodaljac(p,u); % Jac of 'nonlinearity' 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];% put f_u into block matrix 
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
D=u(p.nu+4); Ks=p.mat.K; K=[Ks 0*Ks; 0*Ks D*Ks];% diffusion matrix 
Gu=K-p.mat.M*Fu;                             % Jacobian 