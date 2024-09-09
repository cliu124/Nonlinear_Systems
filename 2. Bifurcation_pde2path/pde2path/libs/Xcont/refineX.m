function p=refineX(p,varargin)
% refineX: mesh-adaption for surface meshes, based on p.fuha.e2rs; 
% controlled by p.sw.rlong  
% p=refineX(p)     or      p=refineX(p,sig)
if size(p.tri,2)==2; p=refineX1D(p,varargin{:}); return; end 
sig=p.nc.sig;  if nargin>1; sig=varargin{1}; end
try; DIR=p.DIR; catch; DIR=[]; end
try; NEU=p.NEU; catch; NEU=[]; end
try; rlong=p.sw.rlong; catch; rlong=0; end
try; nobdref=p.sw.nobdref; catch; nobdref=1; end
ids=p.fuha.e2rs(p,sig); 
if nobdref; ids=rmbdtrifromlist(p.tri,ids); end  % rm bd triangles 
if rlong; [Xn,p.tri,p.DIR,p.NEU]=Trefinelong(p.X,p.tri,DIR,NEU,ids); 
else [Xn,p.tri,p.DIR,p.NEU]=TrefineRGB(p.X,p.tri,DIR,NEU,ids); % refine 
end
fprintf(['refineX, ' func2str(p.fuha.e2rs) ', old,new size(X)=%i,%i\n'],p.np,size(Xn,1)); 
% interpolate old u and tau to new triangulation 
nuo=p.nu; npo=p.np; np=size(Xn,1); nt=size(p.tri,1); par=getaux(p); 
uf=p.mat.fill*p.u(1:nuo); upf=p.mat.fill*p.up(1:nuo); 
tauf=p.mat.fill*p.tau(1:nuo); dtau=p.tau(p.nu+1:end); 
neq=p.nc.neq; un=zeros(np*neq,1); upn=un; taun=un; 
for i=1:p.nc.neq
 upn((i-1)*np+1:i*np)=p3interpol(Xn(:,1),Xn(:,2),Xn(:,3),upf((i-1)*npo+1:i*npo),p.X(:,1),p.X(:,2),p.X(:,3),p);    
 un((i-1)*np+1:i*np)=p3interpol(Xn(:,1),Xn(:,2),Xn(:,3),uf((i-1)*npo+1:i*npo),p.X(:,1),p.X(:,2),p.X(:,3),p); 
 p.sw.ips=2;
 taun((i-1)*np+1:i*np)=p3interpol(Xn(:,1),Xn(:,2),Xn(:,3),tauf((i-1)*npo+1:i*npo),p.X(:,1),p.X(:,2),p.X(:,3),p);
 p.sw.ips=2;
end 
p.up=[upn;par]; % for plotting 
p.u=[un; par]; p.tau=[taun; dtau]; p.nt=nt; p.nu=p.nc.neq*np; p.np=np; p.X=Xn;  
% assume we have at most one boundary, and thats coded in p.idx 
b=boundary_faces(p.tri); ids=unique([b(:,1);b(:,2)]); p.idx=ids; 
if size(p.mat.fill,1)~=1 % update fill and drop 
 q=p; pde.grid.p=p.X'; pde.grid.t=p.tri'; pde.grid.e=b'; q.pdeo=pde; p.sw.orgper=1;  
 p=box2per(q,p.sw.bcper); %p.up=[p.mat.drop*upn; par]; 
 taun=p.mat.drop*taun; p.tau=[taun; dtau];
 if p.sw.Xfill~=0;  % update Xfill matrices 
  p.idx1=find(abs(p.X(:,1)-max(p.X(:,1)))<1e-4 | abs(p.X(:,1)-min(p.X(:,1)))<1e-4); 
  p.idx2=find(abs(p.X(:,2)-max(p.X(:,2)))<1e-4 | abs(p.X(:,2)-min(p.X(:,2)))<1e-4);
  p.idx3=find(abs(p.X(:,3)-max(p.X(:,3)))<1e-4 | abs(p.X(:,3)-min(p.X(:,3)))<1e-4);  
  [p.Xfillx,p.Xfilly,p.Xfillz]=Xfillmat(p);
 end 
end  % if there are several (different) boundaries (e.g., demo libri), identify them outside of refineX 
%fprintf('plotting ...\n'); 
pplot(p,11); 

