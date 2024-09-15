close all; keep pphome; global s0 s1 u0 u1 Psi par xi um1 um2 sig;  
%% Preparations: put filenames into fn, set some bvp parameters 
opt=[]; opt.t1=10; opt.nti=31; flip=0; 
sd0='c1'; sp0='pt0'; sd1='c1'; sp1='pt0'; sname='to-c1-0.mat'; 
%sd0='c1'; sp0='pt0'; sd1='c1'; sp1='pt30'; sname='to-c1-29.mat'; opt.nti=51;
fn=setfnflip(sd0,sp0,sd1,sp1,flip); opt=ocstanopt(opt);
opt.rhoi=1; opt.start=1; opt.tv=[];  opt.retsw=0; 
% the solve and continue call, and some plots 
sol=[]; alvin=[0.025 0.05]; v=[-50,30]; opt.msw=1; opt.AbsTol=1e-3; 
opt.Stats_step='on'; opt.vsw=1; opt.Itnlmax=10; opt.Nmax=200;
ta=tic; [alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,sol,[],opt,fn); toc(ta)
sol0=sol; 
%% cp and diagnostics plot
tit=['Eq to ' fn.sd1 '/' fn.sp1]; pcmp=1; v=[-30,30];
%tit=[fn.sd0 '/' fn.sp0 ' to ' fn.sd1 '/' fn.sp1];
psol3D(s1,sol,4,pcmp,tit,0); view(v); colormap cool; 
switch pcmp; case 1; zlabel('v_1'); case 2; zlabel('v_2'); 
    case 3; zlabel('\lambda_1');  case 4; zlabel('\lambda_2'); end
zdia=lvdiagn(sol,16,fn,1); set(gca,'XTick',[0.1 1 10]);
%% a subsequent call to iscnat
opt.tv=sol.x; opt.start=0; alvin=[0.3 0.6 1]; opt.msw=0; opt.vsw=1;
[alv1,vv1,sol,udat,tlv1,tv1,uv1]=iscnat(alvin,sol,udat.usec,opt,fn); 
alv=[alv alv1]; vv=[vv vv1]; tlv=[tlv tlv1]; tv=[tv; tv1]; uv=[uv; uv1]; 
sol1=sol; 
%% save a path 
save(sname,'fn','Psi','xi','alv','vv','tlv','sol','tv','uv','udat','sig','um1','um2','opt'); 
%% load again (typically for later plotting) 
[fn,alv,vv,tlv,sol,tv,uv,opt]=loadcp(sname);  