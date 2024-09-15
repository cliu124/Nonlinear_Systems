function slsolplot(p,v) 
sol=p.cp; dia=sldiagn(p,15); s1=p.oc.s1; 
psol3DT(s1,sol,1,1,v,[]); zlabel('v'); xlabel('x'); ylabel('t'); % plot P 
psol3DT(s1,sol,2,0,v,[]); zlabel('q'); xlabel('x'); ylabel('t'); % plot k 
