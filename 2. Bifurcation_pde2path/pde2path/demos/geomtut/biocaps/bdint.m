function I=bdint(p,f) % boundary integral, trapz rule
bde=p.DIR; X=p.X; iv=bde(:,1); jv=bde(:,2); 
edl=vecnorm(X(jv,:)-X(iv,:),2,2); 
I=0.5*sum(edl.*(f(iv)+f(jv))); 
