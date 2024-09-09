%% Basic canonical paths of ODE toy system, Cell1, initialization 
p=[]; poc=[]; nt1=40; nt2=4*nt1;  % # of time-slices for CPS and CP (initial) 
om=1; rho=1; th=1; par=[om;rho;th]; % parameters
p=toyinit(p,2,nt1,par); % construction of explicit CPS; now initialize CP data, 
% with mtom (0/1), end-error (derr) and length wT in multiples of CP-length
mtom=1; tadevs=1e-4; poc=ocinit_sp(poc,p,[4;0;0;0],nt2,mtom,tadevs); 
%% Cell 2: CP from states (4,0) to the CPS at states (1,0)
% Here no time adaption is necessary, i.e. this is the easiest case 
poc2=isc(poc,0:0.1:1); toyplot(poc2); 
%% Cell 2b: plot cp in (y1,r), and add some candy for Fig in tutorial 
figure(10); clf; plot(1,1,'b*'); hold on; s=5; oy=0.1; ox=0.03; 
[y1,r]=plotpp(poc2); hold on; z1=1./sqrt(y1); plot(y1,z1,'r'); y1l=length(y1); 
legend('CPS','CP','1/sqrt(y_1)');  xlabel('y_1'); ylabel('r'); 
ar1=annotation('arrow'); ar1.X=0.8*[0.7 0.8]; ar1.Y=0.05*[1 1]; 
plot(y1(1:s:y1l)+ox,z1(1:s:y1l)+oy,'r^','HandleVisibility','off'); 
plot(y1(1:s:30)-ox,z1(1:s:30)-oy,'rv','HandleVisibility','off'); 
axis([0 1 0 10]); set(gca,'fontsize',12); hold off; 
%% Cell 3: CP from states (4,4) to the CPS at states (1,0). The path only 
% needs 7/8 of a circle, i.e. truncation time needs to be adapted. 
poc3=poc; poc3.oc.s0.u(1:4)=[4;4;0;0]; % overwrite starting point
poc3.oc.mtom=0; poc3=isc(poc3,0:0.2:1); toyplot(poc3);
%% Cell 4: CP from states (4,0) to the CPS with different base-point
poc4=poc; s1=poc4.oc.s1; 
s1.hopf.y(1:2,1:end-1)=circshift(s1.hopf.y(1:2,1:end-1),300,2);   % shift CPS
s1.hopf.y(1:2,end)=s1.hopf.y(1:2,1); poc4.oc.s1=s1; 
poc4=isc(poc4,0:0.1:1); toyplot(poc4);