function H=getmeancurv(p,X) % mean curv based on 1st and 2nd fundamental form 
[E,F,G]=get1ff(p,X); [L,M,N]=get2ff(p,X); H=0.25*(L.*G-2*M.*F+N.*E)./(E.*G-F.^2);    