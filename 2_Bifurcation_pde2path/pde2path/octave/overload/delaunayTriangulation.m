function dt=delaunayTriangulation(x)  % emulate matlab's delaunay 
t=delaunayn(x); t=elem_fix(t,x); e=freebdHU(x,t); 
dt.ConnectivityList=t; dt.freeBoundary=e; dt.Points=x; 
end