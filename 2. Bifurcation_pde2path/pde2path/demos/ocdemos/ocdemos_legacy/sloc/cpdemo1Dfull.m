% driver script for SLOC canonical paths 
global s0 s1 u0 u1 Psi par xi um1 um2 sig;  % set global vars 
%% Preparations: put filenames into fn (choose one), set some bvp parameters 
sd0='f1'; sp0='pt13'; sd1='p3'; sp1='pt19'; flip=1; % p3->FSC
sd0='f2'; sp0='pt12'; sd1='p3'; sp1='pt19'; flip=1; % p3->FSM
%sd0='p1'; sp0='pt68'; sd1='p3'; sp1='pt19'; flip=1; % p3->PS
fn=setfnflip(sd0,sp0,sd1,sp1,flip); opt=[]; opt=ocstanopt(opt);
opt.rhoi=1; opt.t1=100; opt.start=1; opt.tv=[]; opt.nti=10; opt.retsw=0; 
% the solve and continue call, and some plots 
sol=[]; alvin=[0.1 0.25]; v=[15,30]; opt.msw=1; opt.Itnlmax=10;
[alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,sol,[],opt,fn); slsolplot(sol,v); 
%% a subsequent call to iscnat
opt.tv=sol.x; opt.start=0; alvin=[0.25 0.5 0.7 1];  opt.msw=0; opt.vsw=0;
[alv1,vv1,sol,udat,tlv1,tv1,uv1]=iscnat(alvin,sol,udat.usec,opt,fn); 
alv=[alv alv1]; vv=[vv vv1]; tlv=[tlv tlv1]; tv=[tv; tv1]; uv=[uv; uv1];
slsolplot(sol,v);
%% values; 
rho=s1.u(s1.nu+1); jp=jcai(s1,sol,rho)+disjca(s1,sol,rho); % jca of path
j0=jca(s0,s0.u)/rho; j1=jca(s1,s1.u)/rho; al=alv(end); % jcas of CSS
fprintf([fn.sd0 '/' fn.sp0 ' to ' fn.sd1 '/' fn.sp1 ', al=' mat2str(al,3) ': ']); 
fprintf('J0=%g, Jp=%g, J1=%g\n',j0,jp,j1); zdia=sldiagn(sol,15);
tit=['J_0=' mat2str(j0,5) ',  J=' mat2str(real(jp),5) ',  J_1=' mat2str(j1,5)]; title(tit); 
%% a simple plot of J over alpha  
figure(6); clf; plot(alv(1,:),vv(1,:),'-*'); set(gca,'FontSize',s1.plot.fs); 
xlabel('\alpha','FontSize',s1.plot.fs); ylabel('J_{a}','FontSize',s1.plot.fs);
%% ---- A fold in alpha, here iscarc needed; Prep. and initial iscarc call
sd0='f1'; sp0='pt13'; sd1='p1'; sp1='pt68'; flip=1; fn=setfnflip(sd0,sp0,sd1,sp1,flip); 
esol=[]; usec=[]; opt.nsteps=3; opt.alvin=[0.2 0.3]; sig=0.1; opt.nti=10; opt.tv=[];
opt.Stats_step='on'; opt.start=1; opt.sigmax=1; opt.retsw=1;
[alv,vv,usec1,esol1,tlv,tv,uv]=iscarc(esol,usec,opt,fn); opt.start=0; 
%% subsequent iccarc-calls (repeat this cell) 
opt.nsteps=35; usec=usec1; esol=esol1; % new input (for repeated calls) 
[alv1,vv1,usec1,esol1,tlv1,tv1,uv1]=iscarc(esol,usec,opt,fn); 
% append output to prev. steps for repeated calls 
alv=[alv alv1]; vv=[vv vv1]; tlv=[tlv tlv1]; tv=[tv; tv1]; uv=[uv; uv1]; 
%% save results for skibademo.m 
alv0=alv; vv0=vv; tv0=tv; uv0=uv; tlv0=tlv; 
%% a simple plot of J over alpha  
figure(6); clf; plot(alv(1,:),vv(1,:),'-*'); set(gca,'FontSize',s1.plot.fs); 
xlabel('\alpha','FontSize',s1.plot.fs); ylabel('J_{a}','FontSize',s1.plot.fs);
%% save a path 
usec=usec1; sols=esol1; % only update (for saving) after success in iccarc! 
save('p68toFSC','fn','Psi','xi','alv','vv','tlv','sols','tv','uv','usec','sig','um1','um2','opt'); 
%% load again (typically for later plotting) 
[fn,alv,vv,tlv,esol1,tv,uv,usec1,opt]=loadcp('p51toFSC');  
%% evaluate selected path from continuation
j=length(tlv); 
tl=tlv(j); v=[25,15]; n=s1.nu; sol.x=tv(j,1:tlv(j));sol.y=squeeze(uv(j,1:n,1:tlv(j)));
al=alv(j), u0=al*s0.u(1:n)+(1-al)*s1.u(1:n); u1=s1.u(1:s1.nu); slsolplot(sol,v); 
%% values; 
rho=s1.u(s1.nu+1); jp=jcai(s1,sol,rho)+disjca(s1,sol,rho); % jca of path
j0=jca(s0,s0.u)/rho; j1=jca(s1,s1.u)/rho; al=alv(end); % jcas of CSS
fprintf([fn.sd0 '/' fn.sp0 ' to ' fn.sd1 '/' fn.sp1 ', al=' mat2str(al,3) ': ']); 
fprintf('J0=%g, Jp=%g, J1=%g\n',j0,jp,j1); zdia=sldiagn(sol,15);
tit=['J_0=' mat2str(j0,5) ',  J=' mat2str(real(jp),5) ',  J_1=' mat2str(j1,5)]; title(tit); 
%% fix al from iscarc to some given value and compute CPs
j=22; alv(j), n=s1.nu; sol.x=tv(j,1:tlv(j));sol.y=squeeze(uv(j,1:n,1:tlv(j))); 
al=0.6; u0=al*s0.u(1:n)+(1-al)*s1.u(1:n); u1=s1.u(1:s1.nu); 
opt.M=s1.mat.M; sol=mtom(@mrhs,@cbcf,sol,opt); v=[100,30]; slsolplot(sol,v); 
%% now the same on lower branch 
j=34; alv(j), pause; n=s1.nu; sol.x=tv(j,1:tlv(j));sol.y=squeeze(uv(j,1:n,1:tlv(j))); 
al=0.6; sol=mtom(@mrhs,@cbcf,sol,opt);  slsolplot(sol,v); 



