function p=oosetfemops(p) % 2D, chebychev-Lap with NBCs 
nx=p.nx; ny=p.ny; p.mat.M=speye(2*nx*ny); g=[0 1 0; 0 1 0]; 
[x,D2x]=cheb2bc(nx,g); [y,D2y]=cheb2bc(ny,g); % weideman 
p.x=x; p.y=y; [xx,yy]=meshgrid(x,y); xx=xx(:); yy=yy(:); 
lb=find(xx==-1); rb=find(xx==1); % left and right boundary 
bb=find(yy==-1); ub=find(yy==1); % bottom and top boundary 
p.bdi=[lb;rb;bb;ub]; % boundary indizes 
L=kron(D2x,eye(ny))./p.lx^2+kron(eye(nx),D2y)./p.ly^2; % Lapl. 
p.mat.L=sparse(L); p.x=x; p.y=y; 