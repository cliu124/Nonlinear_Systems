function p=oosetfemops(p) % 2D, generate Chebychev-Lap with NBCs 
nx=p.nx; ny=p.ny; p.mat.M=speye(nx*ny);
[Dx,x]=chebr(nx+1); [Dy,y]=chebr(ny+1); % Diff.matr.with two add points for BCs 
p.x=x; p.y=y; [xx,yy]=meshgrid(x,y); xx=xx(:); yy=yy(:); %p.xx=xx; 
lb=find(xx==-1); rb=find(xx==1); % left and right boundary 
bb=find(yy==-1); ub=find(yy==1); % bottom and top boundary 
p.bui=setdiff(1:(nx+2)*(ny+2),[lb;rb;bb;ub]); % bulk indizes 
try LL=load(p.Lfn,'LL'); p.mat.L=LL.LL; % try to load L from disk 
catch; % generate L via numjac for the contribution of the boundary modes 
global pj % use a global variable to keep the numjac interface simple 
D2x=Dx^2; D2y=Dy^2; 
L=kron(D2x,eye(ny+2))./p.lx^2+kron(eye(nx+2),D2y)./p.ly^2; % Lapl. on full dom 
p.mat.L=L; p.x=x; p.y=y; pj=p; 
nn=nx*ny; u=zeros(nn,1); r=u; thresh=1e-6*ones(nn,1); 
[Gu,njfac,njG]=numjac('lres',0,u(:),r,thresh,[],0,[],[]); % brute force numjac 
LL=sparse(Gu); % Laplacian only on active nodes 
save(p.Lfn,'LL'); p.mat.L=LL; 
end    % of try-catch 