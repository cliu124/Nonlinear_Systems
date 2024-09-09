function p=oosetfemops(p) % 2D, cheb, with DBCs 
nx=p.nx; ny=p.ny; p.mat.M=speye(nx*ny); [Dx,x]=chebr(nx+1); 
[Dy,y]=chebr(ny+1); % Diff. matrices with two extra points for BCs 
p.x=x; p.y=y; [xx,yy]=meshgrid(x,y); xx=xx(:); yy=yy(:); 
p.lb=find(xx==-1); p.rb=find(xx==1); % left and right Bdry
p.bb=find(yy==-1); p.ub=find(yy==1); % bottom and top Bdry
p.bui=setdiff(1:(nx+2)*(ny+2),[p.lb;p.rb;p.bb;p.ub]); % bulk indizes 
D2x=Dx^2; D2y=Dy^2; 
p.mat.L=kron(D2x,eye(ny+2))/p.lx^2+kron(eye(nx+2),D2y)/p.ly^2; % Laplacian 