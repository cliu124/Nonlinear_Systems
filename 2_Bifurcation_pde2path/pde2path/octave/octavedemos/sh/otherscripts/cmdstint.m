%% init and zero-branch 
lx=4*pi; nx=round(3*lx);ly=2*pi/sqrt(3); lam=-0.001; nu=1.3; par=[lam; nu];  sw.sym=2; sw.ref=1; 
ndim=2; p=shinit(p,nx,lx,ly,ndim,par,sw); p.np, dir='tint'; p=setfn(p,dir); huclean(p); 
po=getpte(p); x=po(1,:)'; y=po(2,:)'; % extract coord from p 
%% a 'good' initial guess for hex2str front, and hex2str front from Newton loop 
p.u(1:p.np)=(cos(x)+cos(x/2).*cos(sqrt(3)*y/2)).*(x<=0)+2*cos(x).*(x>0); p.u(p.nu+1)=0.2; 
spl(p,''); title('initial guess 1'); r=norm(resi(p,p.u),'inf'); 
[u,res,iter]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu); 
fprintf('initial res=%g, res=%g after %i iteration\n',r,res,iter); 
spl(p,''); title('solution 1'); 
%% a 'bad' initial guess for hex2str front, Newton loop goes to stripes 
p.u(1:p.np)=(cos(x)+cos(x/2).*cos(sqrt(3)*y/2)).*(x<=0)+4*cos(x).*(x>0); 
spl(p,''); title('initial guess 2'); r=norm(resi(p,p.u),'inf'); 
[u,res,iter,Gu,Glam,p]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu); 
fprintf('initial res=%g, res=%g after %i iteration\n',r,res,iter); 
spl(p,''); title('solution 2'); 
%% a 'bad' initial guess for hex2str front, do a few steps with tintxs
p.u(1:p.np)=(cos(x)+cos(x/2).*cos(sqrt(3)*y/2)).*(x<=0)+4*cos(x).*(x>0); 
t1=0; ts=[]; nc=0; dt=0.01; nt=500; pmod=10; smod=100; p.mat.Kadv=0; 
%% the tint loop, repeat this cell until residual is small (here just once) 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@nodalf); 
%% Newton loop after tint 
[u,res]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu); res, plotsol(p,1,1,2); 
%% plot time series of res 
tss=5; plot(ts(1,tss:end),ts(2,tss:end)); axis tight; legend('res'); xlabel('t'); 
%% cont of the obtained solution
clf(2);p=pmcont(p,20);
