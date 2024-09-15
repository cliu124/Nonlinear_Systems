format compact; close all; keep pphome; p=[];
%%
close all; p=stanparam(p); screenlayout(p);
p.nc.dsmin=1e-6; p.nc.dsmax=0.1; p.nc.imax=20;
p.nc.tol=1e-10; p.nc.nsteps=10; p.nc.neig=10; 
p.nc.neq=2; p.nc.ilam=1; p.nc.nq=0; p.sw.spcalc=0; p.sw.bprint=[]; 
p.sw.jac=1; p.sw.para=1; p.nc.lamdtol=1e-1;
p.plot.pcmp=2; p.plot.bpcmp=2; p.plot.cm='cool';p.plot.pstyle=3;
p.nc.lammin=-3; p.nc.lammax=6;% p.nc.del=1e-3; 
p.file.smod=5; lx=pi; nx=40; 
p.nc.dlammax=0.1; p.sw.bifcheck=0; p.pffac=0;
kstar=[0.5 0.5]; % wave vector of the solution
p.sig=-1;om=0.46; % PDE coefficient sigma and a selected appr. omega
p.sw.bcper=[1 2];   % periodicity: in both the x and y directions
icsw=0; % initial guess switch: 0 for u=0, 1 for asymptotics based on the eifenfunction
p=nlbinit(p,lx,nx,om,kstar,icsw); 
%% detect bifurcation points along the trivial branch, using findbifm 
% (i.e. findbif modified to detect double evals crossing the imag. axis)
p=setlam(p,0.16); % set lam to startvalue for bif. point search
p=resetc(p); p.u(1:p.nu)=zeros(p.nu,1);
p.nc.bisecmax=15; p.nc.nsteps=30;p.nc.dsmin=1e-12;
p.sol.ds=0.05;p.sol.restart=1; clf(2)
p=findbifm(p,4); % detect bif. point(s)
%% branch switching, if pffac (phase-fix-factor) >0, then fix phase 
nr=1; % use first e-vec of the two at bif point (unless user choses otherwise)
% resulting branch is independent of nr as it should be (for pffac=0) 
% For phase-fix set pffac<0 for dropp, or pffac>0 for "method (b)" 
pffac=1e3; del=1e-12; 
q=swibram('p','bpt1','q1',0.02,del,nr,pffac);
q.nc.nsteps=10; q=cont(q);
q=swibram('p','bpt2','q2',0.02,1e-12,nr,pffac);
q.nc.nsteps=10; q=cont(q);
q=swibram('p','bpt3','q3',0.02,1e-12,nr,pffac);
q.nc.nsteps=10; q=cont(q);
%% plot BD
figure(3); clf; cmp=2; lms=5; 
plotbraf('q1','pt10',3,cmp,'cl','b','lab',5,'lms',lms,'lwst',2); 
plotbraf('q2','pt10',3,cmp,'cl','r','lab',5,'lms',lms); 
plotbraf('q3','pt10',3,cmp,'cl','m','lab',5,'lms',lms); 
axis([0 0.6 0 2]);
xlabel('\omega'); ylabel('||u||_2');
%% solution plots
a=15;b=65;
plotsolf('q1','pt5',1,1,3);xlabel(''); ylabel('');zlabel('');
set(gca,'ztick',[-0.2 0 0.2]);view(a,b);pause 
plotsolf('q2','pt5',1,1,3);xlabel(''); ylabel('');zlabel('');
set(gca,'ztick',[-0.2 0 0.2]);view(a,b);pause 
plotsolf('q3','pt5',1,1,3);xlabel(''); ylabel('');zlabel('');
set(gca,'ztick',[-0.2 0 0.2]);view(a,b);