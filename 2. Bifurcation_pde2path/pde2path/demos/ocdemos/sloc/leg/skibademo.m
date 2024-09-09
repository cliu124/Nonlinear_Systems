% Skiba example, continues cpdemo1D.m, legacy setup 
% find paths from the uv0 initial states from cpdemo.m to FSM 
js=10; je=30; jl=je-js+1; % alpha-range; now set the target and Psi to FSM: 
s1=loadp('f2','pt12'); u1=s1.u(1:s1.nu); [Psi,muv,d,t1]=getPsi(s1); 
a0l=length(alv0); tva=zeros(jl,opt.Nmax+1); % some prep. and fields to hold paths 
n=s1.nu; uva=zeros(jl,n+1,opt.Nmax+1);  alva=[]; vva=[]; tavl=[]; sol=[]; 
alvin=[0.1 0.25 0.5 0.75 1]; % we run from uv0(j,:) to FSM with iscnat
tv=linspace(0,opt.t1,opt.nti); se=2; opt.tv=tv.^se./opt.t1^(se-1); doplot=1; 
opt.start=0; opt.msw=0; opt.Stats_step='off'; v=[50,8]; % switch off stats 
for j=js:je; % continue till Skiba paths found, or till the end ...
   fprintf('j=%i, al=%g\n', j, alv0(j)); s0.u(1:n)=uv0(j,1:n,1)'; 
   [alv,vv,sol,udat]=iscnat(alvin,[],[],opt,fn); 
   if alv(end)==1; Jd=vv0(j)-vv(end); fprintf('J1-J2=%g\n',Jd); % cont. success 
      alva=[alva alv0(j)]; vva=[vva vv(end)]; tl=length(sol.x); % save values 
      talv=[tavl tl]; tva(j,1:tl)=sol.x; uva(j,1:n,1:tl)=sol.y; 
      if abs(Jd)<0.05; doplot=asknu('plot path?',doplot); % Skiba point(s) found
        if doplot==1; sol0=[]; alp=alv0(j); % plot the paths to FSC and FSM 
          sol0.x=tv0(j,1:tlv0(j));sol0.y=squeeze(uv0(j,1:n,1:tlv0(j))); 
          psol3Dm(s1,sol0,sol,1,1,[]); view(v); zlabel('P'); % plot both paths
          psol3Dm(s1,sol0,sol,2,0,[]); view(v); zlabel('k'); pause
        end
      end
   end
end
%% plot value diagram 
jep=je; figure(6); clf; plot(alv0(js:jep),vv0(js:jep),'-*b');hold on;
plot(alva,vva,'-*r');set(gca,'FontSize',14); axis tight;
xlabel('\alpha','FontSize',14); ylabel('J_{a}','FontSize',14); 