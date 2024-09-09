function Gu=sGjac1D(p,u)  % Jac for ql-AC 
par=u(p.nu+1:end); lam=par(2); ga=par(3); del=par(4); epsi=par(5); 
n=p.nu; u=u(1:n); M=p.mat.M; gr=p.pdeo.grid; fem=p.pdeo.fem;  
fu=lam+3*u.^2-5*ga*u.^4; Fu=spdiags(fu,0,n,n); % Fu as usual 
ut=p.mat.p2c*u; c=cfu(ut,par);  [K,~,~]=fem.assema(gr,c,0,0);
switch p.jacsw; 
    case 1; % approximate version
        cu=del+2*epsi*u; ux=p.mat.Dx*u; % 1st derivative as coefficient
        K1=p.mat.Kx*spdiags(cu.*ux,0,n,n); % first order derivative acting on v 
        Gu=K-K1-M*Fu+p.nc.sf*p.mat.Q; % putting it all together
    case 2;  % numjac for \pa_u div(c(u)\nab v)    
        Kuvd=getKuvd(p,par,u,u); Gu=K+Kuvd-M*Fu+p.nc.sf*p.mat.Q; 
end