function mypsol(dir,pt); 
p=loadpp(dir,pt); p.plot.shsw=1; p.plot.cm='hot'; plotsol(p,1,1,3); 
xticks([-2 2]); yticks([-2 2]); nola; zlabel(''); 