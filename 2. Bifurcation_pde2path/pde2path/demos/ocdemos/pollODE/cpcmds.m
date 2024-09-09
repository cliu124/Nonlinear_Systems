%% ODE-pollution CPs; 
% first CP from [0 0] to the CSS at rho=0.55, T-adaptation works nicely 
fdir='FSS'; fpt='pt7'; 
p1=[]; p1=ocinit(p1,fdir,fpt,fdir,fpt); % load standard options
p1.oc.nti=500; p1.oc.rhoi=1;  p1.oc.mtom=0; p1.oc.tadevs=1e-1; 
p1.tomopt.maxIt=20; p1.oc.T=800; 
p1.u0(1:2)=[0 0]; alv1=[0.01 0.05 0.1]; alv2=0.2:0.2:1; 
tic; p2a=isc(p1,alv1); toc, cpplot(p2a,4); % go and plot 
%% go further 
p2b=isc(p2a,alv2); cpplot(p2b,4); tadev(p2b); 
%% CP to CSS with u0=[0.4 0.4], again 2 steps 
p1.u0(1:2)=[0.4 0.4]; alv1=[0.05 0.1]; p3=isc(p1,alv1); cpplot(p3,4); pause; 
alv2=0.2:0.1:1; p4=isc(p3,alv2); cpplot(p4,4); 
%% CP to the hom.CPS at rho=0.57; Preparations: load CPS and optionally refine 
hdir='h1'; hpt='pt13'; fdir='FSS'; fpt='pt11'; 
[muv1, muv2,ind,h]=floqpsap(hdir,hpt); fprintf('d(u_H)=%i\n',ind); 
q=loadp(hdir,hpt); huclean(q); hoplot(q,1,1);  % load hopf orbit, rho=0.57
if 0 % optionally refine hopf orbit in time; here not needed 
q=uhopftref(q,1.5); % refine mesh
[y,T,lam,res,iter,A,q]=honloop(q,q.hopf.y,q.hopf.T,q.u(q.nu+q.nc.ilam)); 
q.hopf.y=y; q.hopf.T=T; res, hoplot(q,2,1); end
%% CP from [v,w]=[0.4, 0.4] to hom CPS, init 
p1=[]; p1=ocinit(p1,fdir,fpt,hdir,hpt); % load standard options
p1.oc.s1=q; p1.oc.rhoi=1; p1.oc.mtom=0; 
p1.oc.tadevs=1e-3; p1.oc.nti=100; % start with small nti
p1.u0(1:2)=[0.4 0.4]; al1=0.1; al2=0.2:0.2:1; 
p2a=isc(p1,al1); cpplot(p2a,4); % go and plot 
%% go further, 
p2b=isc(p2a,al2); cpplot(p2b,4); tadev(p2b); 
%% CP from last initial states to shifted hom CPS, first shift 
q=loadp(hdir,hpt); y=q.hopf.y; tl=q.hopf.tl; 
q.hopf.y=circshift(y(:,1:end),round(tl/2),2); 
[y,T,lam,res,iter,A,q]=honloop(q,q.hopf.y,q.hopf.T,q.u(q.nu+q.nc.ilam)); 
q.hopf.y=y; q.hopf.T=T; res, hoplot(q,2,1); p1.oc.s1=q; 
%% compute initial CP 
p4a=isc(p1,0.25); cpplot(p4a,4); p4a.oc.tadevs=1e-2; 
%% go further 
p4b=isc(p4a,0.4:0.2:1); cpplot(p4b,4); tadev(p4b); 
%% a-posteriori decrease of target deviation 
fdir='FSS'; fpt='pt2'; p=[]; p1=ocinit(p1,fdir,fpt,fdir,fpt);  
p1.oc.nti=1000; p1.oc.rhoi=1;  p1.oc.mtom=0; p1.oc.tadevs=0.1; p1.oc.T=200; 
p1.u0(1:2)=[0 0]; alv1=[0.01 0.05 0.1]; alv2=0.2:0.2:1; 
p2a=isc(p1,alv1); cpplot(p2a,4); pause; % initial steps 
p2b=isc(p2a,alv2); cpplot(p2b,4); tadev(p2b); % ini
%% decrease ||u(1)-uh0||
p2c=p2b; p2c.oc.tadev2=5e-3; p2c=isc(p2c,1); cpplot(p2c,4); tadev(p2c); 