function c=cfu(u,par) % diffusion tensor, for \pa_u div(c(u)grad v)
c0=par(1); del=par(4); epsi=par(5); 
c=c0+del*u+epsi*u.^2; 
end 