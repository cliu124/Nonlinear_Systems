function djca=disjcaT(p,sol,rho,s1ho) 
% disjcaT: discounted current value jca 
par=p.u(p.nu+1:end); tl=length(sol.t); t=sol.par(1);
if s1ho==0 % CSS / CP to CSS
    djca=jca(p,[sol.u(:,tl);par])*exp(-rho*t)/rho;
else % CPS / CP to CPS
    T=p.hopf.T; sol2.t=p.hopf.t;  sol2.u=p.hopf.y;
    sol2.par=p.hopf.T; C=jcaiT(p,sol2,rho);
    djca=C*exp(-rho*t)/(1-exp(-T));
end

