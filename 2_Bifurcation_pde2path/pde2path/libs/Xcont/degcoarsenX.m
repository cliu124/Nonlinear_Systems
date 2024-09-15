function p=degcoarsenX(p,sig,varargin) 
% degcoarsenX:  coarsen via rmdegfaces 
%
% p=degcoarsenX(p,sig);           % sig=fraction of triangles to coarsen 
% p=degcoarsenX(p,sig,it);        % it=#iterations (default 5)
% p=degcoarsenX(p,sig,it,keepbd); % keepd boundary triangles for keepbd=1 (default 0)
% 
it=5; keepbd=0; if nargin>2; it=varargin{1}; end 
if nargin>3; keepbd=varargin{2}; end 
try areafac=p.nc.areafac; catch; areafac=1; end
np=p.np; ic=max(round(sig*np),1); % ic=10; 
A=doublearea(p.X,p.tri); [As,idx]=sort(A,'ascend'); c=1; % internal coarsening control 
epsi=c*As(ic); [Xn,p.tri]=rmdegfaces(p.X,p.tri,'MaxIter',it,'Epsilon',epsi,'keepbd',keepbd,'areafac',areafac); 
fprintf('degcoarsenX, old,new size(X)=%i,%i\n',np,size(Xn,1)); 
nuo=p.nu; npo=p.np; np=size(Xn,1); p.nt=size(p.tri,1); par=getaux(p);  
uf=p.mat.fill*p.u(1:nuo); upf=p.mat.fill*p.up(1:nuo); 
tauf=p.mat.fill*p.tau(1:nuo); dtau=p.tau(nuo+1:end); 
neq=p.nc.neq; un=zeros(np*neq,1); upn=un; taun=un; % interpol sol/tangent to new mesh 
for i=1:p.nc.neq
 un((i-1)*np+1:i*np)=p3interpol(Xn(:,1),Xn(:,2),Xn(:,3),uf((i-1)*npo+1:i*npo),p.X(:,1),p.X(:,2),p.X(:,3),p); 
 upn((i-1)*np+1:i*np)=p3interpol(Xn(:,1),Xn(:,2),Xn(:,3),upf((i-1)*npo+1:i*npo),p.X(:,1),p.X(:,2),p.X(:,3),p); 
 taun((i-1)*np+1:i*np)=p3interpol(Xn(:,1),Xn(:,2),Xn(:,3),tauf((i-1)*npo+1:i*npo),p.X(:,1),p.X(:,2),p.X(:,3),p); 
end 
p.u=[un; par]; p.tau=[taun; dtau]; p.up=[upn;par]; % for plotting 
p.nu=p.nc.neq*np; p.np=np; p.X=Xn; 
b=boundary_faces(p.tri); ids=unique([b(:,1);b(:,2)]); p.idx=ids; 
% possibly update drop and fill 
if size(p.mat.drop,1)~=1 
 q=p; pde.grid.p=p.X'; pde.grid.t=p.tri'; pde.grid.e=b'; q.pdeo=pde; p=box2per(q);
 %p.up=p.mat.drop*p.up(1:p.np); p.up=[upn;par]; 
 taun=p.mat.drop*taun; p.tau=[taun; dtau];
 if p.sw.Xfill~=0;  % update Xfill matrices 
  p.idx1=find(abs(p.X(:,1)-max(p.X(:,1)))<1e-4 | abs(p.X(:,1)-min(p.X(:,1)))<1e-4); 
  p.idx2=find(abs(p.X(:,2)-max(p.X(:,2)))<1e-4 | abs(p.X(:,2)-min(p.X(:,2)))<1e-4);
  p.idx3=find(abs(p.X(:,3)-max(p.X(:,3)))<1e-4 | abs(p.X(:,3)-min(p.X(:,3)))<1e-4);  
  [p.Xfillx,p.Xfilly,p.Xfillz]=Xfillmat(p);
 end 
end
pplot(p,11); colorbar; 
