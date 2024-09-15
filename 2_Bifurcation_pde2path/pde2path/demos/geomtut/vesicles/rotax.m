function p=rotax(p,m) % find rotational axis for axisym.vesicle 
nX=vecnorm(p.X,2,2); 
if isequal(m,'max'); i=find(nX(1:p.np)==max(nX(1:p.np))); 
elseif isequal(m,'min'); i=find(nX(1:p.np)==min(nX(1:p.np))); end
p.w=normalizerow(p.X(i(1),:)); 
p.om=0.5*[-p.w(2) p.w(1) 0]+0.5*[0 -p.w(3) p.w(2)]+0.5*[-p.w(3) 0 p.w(1)]; 
p.om=normalizerow(p.om); p.rh=cross(p.om,p.w); p.rh=normalizerow(p.rh); 