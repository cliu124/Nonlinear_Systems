%% Skiba-demo for vegOC; 
global s0 s1 u0 u1 Psi par xi um1 um2 sig;  
%% Skiba at R=28: first connect intermediate p1 to FSS 
sd0='FSS'; sp0='pt14'; sd1='p1'; sp1='pt16'; flip=1; % p1i to FSS 
fn=setfnflip(sd0,sp0,sd1,sp1,flip); 
alvin=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; opt.Nmax=100;
opt.start=1; opt.tv=[]; opt.nti=60; opt.retsw=1;  opt.t1=100; 
[alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,[],[],opt,fn); 
%% store
tlv0=tlv; tv0=tv; uv0=uv; alv0=alv; vv0=vv;
%% now connect IC from uv0 to p1/pt34 (upper p1) 
js=8; je=10; jl=js-je+1; % alpha-range; now set the target and Psi: 
s1=loadp('p1','pt34'); getaux(s1)', fn.sd1='p1'; fn.sp1='pt34'; 
u1=s1.u(1:s1.nu); [Psi,muv,d,t1]=getPsi(s1); n=s1.nu;
a0l=length(alv0); tva=zeros(jl,opt.Nmax+1); % some prep. and fields to hold paths 
opt.Nmax=200; opt.Itnlmax=10; opt.nti=30;
uva=zeros(jl,n+1,opt.Nmax+1); opt.start=0;  sol=[]; alva=[]; vva=[]; tavl=[]; 
opt.retsw=0; alvin=[0.05 0.2 0.4 0.6 0.8 1]; 
tv=linspace(0,opt.t1,opt.nti); se=2; opt.tv=tv.^se./opt.t1^(se-1); doplot=1; 
opt.msw=0; opt.Stats_step='on'; v=[10,20];  % switch off stats 
for j=js:je; 
   fprintf('j=%i, al=%g\n', j, alv0(j)); s0.u(1:n)=uv0(j,1:n,1)'; %norm(s0.u(1:n)), pause
   [alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,[],[],opt,fn); %vv(:), pause
   if ~isempty(alv); Jd=vv0(j)-vv(end); fprintf('J1-J2=%g\n',Jd); % cont. success 
      alva=[alva alv0(j)]; vva=[vva vv(end)]; tl=length(sol.x); % put vals in vector
      talv=[tavl tl]; tva(j,1:tl)=sol.x; uva(j,1:n,1:tl)=sol.y; 
   end
end
%% plot value diagram 
js=5; jep=10; 
figure(6); clf; plot(alv0(js:jep),vv0(js:jep),'-*b');hold on;
plot(alva,vva,'-*r');set(gca,'FontSize',s1.plot.fs); axis tight;
xlabel('\alpha','FontSize',s1.plot.fs); ylabel('J_{a}','FontSize',s1.plot.fs);
%% select j and compute path 
alvin=[0.05 0.2 0.4 0.6 0.8 1]
j=8; fprintf('j=%i, al=%g\n', j, alv0(j)); s0.u(1:n)=uv0(j,1:n,1)'; 
[alv,vv,sol,udat,tlv,tv,uv]=iscnat(alvin,[],[],opt,fn); 
%% plot both paths (to FSS and 'upper' p1) into one figure (componentwise) 
sol0=[]; alp=alv0(j); 
sol0.x=tv0(j,1:tlv0(j));sol0.y=squeeze(uv0(j,1:n,1:tlv0(j))); 
psol3Db2(s1,sol0,sol,1,1,[]); xlabel('x'); ylabel('t'); zlabel('v'); view(v); 
psol3Db2(s1,sol0,sol,2,2,[]); xlabel('x'); ylabel('t'); zlabel('w'); view(v); 
psol3Db2(s1,sol0,sol,4,5,[]); xlabel('x'); ylabel('t'); zlabel('E'); view(v);