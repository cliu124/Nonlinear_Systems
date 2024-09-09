function geo=hexgeo(lx,ly) 
% rectangle with upper edge made up of three pieces 
% with middle from -lx to lx and offset ly
if (lx<0)||(lx>=1) fprintf('recx error: lx out of bounds'); end
if (ly<-1) fprintf('recx error: ly out of bounds'); end
x=[-1 -1 1 1 lx -lx]; y=[1 -1 -1 1 1+ly 1+ly]; geo=polygong(x,y);
return 
