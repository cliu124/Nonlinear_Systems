function p=coarsenX(p,varargin) 
% coarsenX: mesh-adaption for surface meshes, based on p.fuha.e2rs; 
%
% p=coarsenX(p,sig);  controled by p.sw.rlong; generally only coarsens 
% triangles that were previously refined; based on Schmidt & Funken 
if size(p.tri,2)==2; p=coarsenX1D(p,varargin{:}); return; end % X=curve
try; DIR=p.DIR; catch; DIR=[]; end
try; NEU=p.NEU; catch; NEU=[]; end
try; nobdcoarsen=p.sw.nobdcoarsen; catch; nobdcoarsen=0; end
try n0=p.n0; catch; n0=200; end 
sig=p.nc.sig;  if nargin>1; sig=varargin{1}; end
ids=p.fuha.e2rs(p,sig);  
if nobdcoarsen; ids=rmbdtrifromlist(p.tri,ids); end  % rm bd triangles 
[Xn,p.tri,p.DIR,p.NEU]=TcoarsenRGB(n0,p.X,p.tri,DIR,NEU,ids); % refine long 
fprintf('coarsenX, old,new size(Xn)=%i,%i\n',p.np,size(Xn,1)); 
% interpolate old u and tau to new triangulation 
nuo=p.nu; npo=p.np; np=size(Xn,1); nt=size(p.tri,1); par=getaux(p);  
uf=p.mat.fill*p.u(1:nuo); upf=p.mat.fill*p.up(1:nuo); 
tauf=p.mat.fill*p.tau(1:nuo); dtau=p.tau(p.nu+1:end); 
neq=p.nc.neq; un=zeros(np*neq,1); upn=un; taun=zeros(np*neq,1); 
for i=1:p.nc.neq
 un((i-1)*np+1:i*np)=p3interpol(Xn(:,1),Xn(:,2),Xn(:,3),uf((i-1)*npo+1:i*npo),p.X(:,1),p.X(:,2),p.X(:,3),p); 
 upn((i-1)*np+1:i*np)=p3interpol(Xn(:,1),Xn(:,2),Xn(:,3),upf((i-1)*npo+1:i*npo),p.X(:,1),p.X(:,2),p.X(:,3),p); 
 taun((i-1)*np+1:i*np)=p3interpol(Xn(:,1),Xn(:,2),Xn(:,3),tauf((i-1)*npo+1:i*npo),p.X(:,1),p.X(:,2),p.X(:,3),p); 
end 
p.u=[un; par]; p.tau=[taun; dtau]; p.up=[upn;par]; % for plotting 
p.nt=nt; p.nu=p.nc.neq*np; p.np=np; p.X=Xn;  
%p=box2per(p); p.up=[p.mat.drop*upn; par];
pplot(p,11); colorbar;
return
if ~isempty(p.DIR) % assume we have at most one boundary, and thats coded in p.idx 
 p.idx=unique([p.DIR(:,1); p.DIR(:,2)]); 
else; q=boundary_faces(p.tri); ids=unique([q(:,1);q(:,2)]); p.idx=ids; 
end  % if there are several (different) boundaries (e.g., demo libri), identify them outside of meshadaX 