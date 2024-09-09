function bc=nbc(p,u) % BC for model rot
b1=10; b2=0.01; % q-vals in \pa_n u+q*u=0 formulation 
c1=p.eqn.c(1); c2=p.eqn.c(end);  % diff. constants 
b1=b1*c1; b2=b2*c2; % org q-vals need to be multipl. by diff-const. 
enum=max(p.mesh.e(5,:)); % #ofboundary segments 
g=[0;0];q=[[b1 0]; [0 b2]];  
bc=gnbc(p.nc.neq,enum,q,g); % same BC on all bdry segments 