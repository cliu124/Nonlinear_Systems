%% CPs for pollution; first to the inhom.CPS; load CPS
hdir='h1'; hpt='pt8'; fdir='FSS'; fpt='pt9'; 
[muv1, muv2,ind,h]=floqpsap(hdir,hpt); fprintf('d(u_H)=%i\n',ind); 
q=loadp(hdir,hpt); huclean(q); hoplot(q,1,1);  % load hopf orbit
if 0 % optionally refine hopf orbit in time; here not needed 
q=uhopftref(q,1.5); % refine mesh
[y,T,lam,res,iter,A,q]=honloop(q,q.hopf.y,q.hopf.T,q.u(q.nu+q.nc.ilam)); 
q.hopf.y=y; q.hopf.T=T; res, hoplot(q,2,1); end
%% init, and initial CP  from [0.205 0.72] to CPS; here ga2=0.9, but conv. good on 
% reasonable t-interval 
p1=[]; p1=ocinit(p1,fdir,fpt,hdir,hpt); p1.oc.s1=q; p1.oc.nTp=20; p1.oc.nti=1500;
p1.oc.mtom=0; ov=ones(q.np,1); p1.u0(1:2*q.np)=[0.205*ov;0.72*ov]; % set initial states 
p2b=isc(p1,0.1); cpplot(p2b,4); % go and plot 
%% go further, first till al=0.75, then further! 
p3b=isc(p2b,0.25:0.25:0.75); cpplot(p3b,4); 
%%
p4b=isc(p3b,0.85:0.05:0.9); cpplot(p4b,4); 
p4c=isc(p4b,[0.95 0.975]); cpplot(p4c,4); 
%%
polldiagn(p4c,11,20,10,1,1); tadev(p4c)
%% CP from states [0, 0] to hom CPS (as in pollODE), with target-control! 
hdir='h2'; hpt='pt17'; fdir='FSS'; fpt='pt13'; 
[muv1, muv2,ind,h]=floqpsap(hdir,hpt); fprintf('d(u_H)=%i\n',ind); 
q=loadp(hdir,hpt);  hoplot(q,1,1); % load hopf orbit
p=[]; p=ocinit(p,fdir,fpt,hdir,hpt); % load standard options
p.oc.s1=q; p.oc.nti=500; p.oc.nTp=2; p.oc.mtom=0; np=q.np; 
ov=ones(np,1); p.u0(1:2*np)=[0*ov;0*ov]; p.oc.tadevs=1e-2;  
p.tomopt.tol=1e-6; p5=isc(p,0.25); cpplot(p5,4); tadev(p5); %  initial CP 
%% go further 
p6=isc(p5,0.5:0.25:1); cpplot(p6,4); polldiagn(p6,11,20,10,1,1); tadev(p6);