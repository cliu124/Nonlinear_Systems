function slsolplot(sol,v) % just to make the script shorter
global s1; udia=sldiagn(sol,15); 
psol3D(s1,sol,1,1,v,[]); zlabel('P'); xlabel('x'); % plot P 
psol3D(s1,sol,2,0,v,[]); zlabel('k'); xlabel('x');% plot k 
