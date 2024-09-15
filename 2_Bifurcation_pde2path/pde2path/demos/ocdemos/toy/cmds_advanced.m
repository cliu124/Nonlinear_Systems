%% Advanced problems for toy-ODE. 
% Cell 1: leading floquet-multipl=0.52, giving slow convergence.
% Small initial T guess leading to T adaption by adding of periods. 
% Optional (projsw=1) analytical projections 
om=4e-2; rho=1; th=1; par=[om;rho;th];  
p=[]; poc=[]; nt1=20; tadevs=1e-3; mtom=0; nt2=60; 
p=toyinit(p,2,nt1,par); poc=ocinit_sp(poc,p,[4;0;0;0],nt2,mtom,tadevs); 
[mon,muv,F,F2]=anaproj(p,p.hopf.y(:,end-1));  projsw=1; muvs=sort(muv); 
fprintf(['analytical multipl:' repmat('%g ',1,4) '\n'], muvs); 
poc.oc.retsw=1; % return history, e.g., for T plot in Cell 2b
if projsw>0; poc.oc.F=F; poc.oc.F2=F2; poc.oc.muv=muv; end             
poc=isc(poc,0.2:0.1:1); toyplot(poc); % continuation call and plotting 
if projsw==0; % show difference between analytical and numerical multipliers 
  fprintf(['analytical multipl:' repmat('%g ',1,4) '\n'], muvs);   
  fprintf(['numerical  multipl:' repmat('%g ',1,4) '\n'], sort(poc.oc.muv));   
end
%% Cell 2: plot deviation from CPS at maxima
[PKS,LOCS]=findpeaks(poc.cp.u(1,:)); maxi=max(poc.oc.s1.hopf.y(1,:));
peakerr=abs(maxi-PKS); npeaks=length(PKS); nplot=12; os=0; 
esterrPeak=peakerr(end-nplot)*poc.oc.muv(3).^(-1*(npeaks:-1:0)+nplot);
figure(5);clf;plot(npeaks-nplot+1:npeaks-os,peakerr(end-nplot+1:end-os));hold on;
plot(npeaks-nplot+1:npeaks-os,esterrPeak(end-nplot+1:end-os));
xlabel('x_1 peak number '); ylabel('Deviation from CPS');
legend('numerical','linear estimate'); 
set(gca,'fontsize',12); axis tight;
%% Cell 2b: T-plot
al=poc.hist.alpha; Tv=0*al; for i=1:length(al); Tv(i)=poc.hist.par{i}; end 
figure(10); plot(al,Tv,'-*'); set(gca,'fontsize',14); axis tight; 
xlabel('\alpha'); ylabel('T'); 