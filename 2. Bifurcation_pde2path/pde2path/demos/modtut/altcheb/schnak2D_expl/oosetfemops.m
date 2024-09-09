function p=oosetfemops(p) % 2D, chebychev-Lap with NBCs 
% here saving to p.LL (value speed over diskspace) 
nx=p.nx; ny=p.ny; p.mat.M=speye(2*nx*ny);
[Dx,x]=chebr(nx+1); [Dy,y]=chebr(ny+1); % Diff. matrices with two extra points for BCs 
p.x=x; p.y=y; [xx,yy]=meshgrid(x,y); xx=xx(:); yy=yy(:); 
lb=find(xx==-1); rb=find(xx==1); bb=find(yy==-1); ub=find(yy==1); % left, right, bottom and top Bdry
p.bui=setdiff(1:(nx+2)*(ny+2),[lb;rb;bb;ub]); % bulk indizes 
%if ~isfield(p,'LL') % generate LL 
try LL=load(p.Lfn,'LL'); p.mat.L=LL.LL; % try to load L from disk 
catch; % generate and save L 
global pj
D2x=Dx^2; D2y=Dy^2; 
L=kron(D2x,eye(ny+2))./p.lx^2+kron(eye(nx+2),D2y)./p.ly^2; % Laplacian on full dom 
% modify L to take NBCs into account
p.mat.L=L;  pj=p; 
nn=nx*ny; u=zeros(nn,1); r=u; thresh=1e-6*ones(nn,1); 
tic; [Gu,njfac,njG]=numjac('lres',0,u(:),r,thresh,[],0,[],[]);toc  % brute force numjac 
LL=sparse(Gu); % Laplacian only on active nodes 
save(p.Lfn,'LL'); 
p.mat.L=LL; 
end 