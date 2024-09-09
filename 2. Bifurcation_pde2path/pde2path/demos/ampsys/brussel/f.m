function w=f(u,p) % nonlinearity (everything except diffusion) for the Brusselator
% u is the vector of unknowns, everything else (e.g., parameters) passed via p
a=p.par(1); b=p.par(2);  % extract parameters
f1=a-(b+1)*u(1)+u(1)^2*u(2); f2=b*u(1)-u(1)^2*u(2);  w=[f1;f2];