close all; clear all; global s0 s1 u0 u1 Psi par xi um1 um2 sig; % set global vars 
%%  driver script for slOC2D  -- this might be expensive
sd0='h2'; sp0='pt16'; sd1='2DFSC'; sp1='pt13'; flip=0; % Fig. 2a
%sd0='h2'; sp0='pt14'; sd1='2DFSC'; sp1='pt10'; flip=0; % OO Test
if flip==1; dud=sd0; dup=sp0; sd0=sd1; sp0=sp1; sd1=dud; sp1=dup; end % FLIP 
fn.sd0=sd0; fn.sp0=sp0; fn.sd1=sd1; fn.sp1=sp1; 
%% choose solver and bvp parameters 
opt=tomset('Stats','off','Stats_step','on', 'Monitor',3,'order',2,'NMax',100); 
opt=tomset(opt,'AbsTol',1e-3,'FJacobian',@fjac, 'Itnlmax',5, 'Itlinmax',10); 
opt=tomset(opt,'BCJacobian',@cbcjac); opt.lu=0; opt.vsw=0; opt.msw=0; 
opt.rhoi=1; opt.t1=50; opt.start=1; opt.tv=[]; opt.nti=10; opt.retsw=0; 
%% the solve and continue call 
sol=[]; alvin=[0.1 0.2]; % first check how long it takes for 2 steps 
[alv,vv,sol,udat,tlv,tv,uv1]=iscnat(alvin,sol,[],opt,fn); 
sol0=sol; zdia=sldiagn(sol,15); % backup sol, and show soln behaviour  
%% a subsequent call to iscnat
opt.tv=sol.x; opt.start=0; alvin=[0.3 0.5 0.75 1]; opt.Stats_step='on'; 
ti1=tic; [alv,vv,sol,udat,tlv,tv,uv1]=iscnat(alvin,sol,[],opt,fn); 
ti1=toc(ti1)
%% evaluate path 
rho=s1.u(s1.nu+1); jp=jcai(s1,sol,rho)+disjca(s1,sol,rho); % jca of path
j0=jca(s0,s0.u)/rho; j1=jca(s1,s1.u)/rho; al=alv(end); % jcas of CSS
fprintf([fn.sd0 '/' fn.sp0 ' to ' fn.sd1 '/' fn.sp1 ', al=' mat2str(al,3) ': ']); 
fprintf('J0=%g, Jp=%g, J1=%g\n',j0,jp,j1); zdia=sldiagn(sol,15);
tit=['J_0=' mat2str(j0,5) ',  J=' mat2str(real(jp),5) ',  J_1=' mat2str(j1,5)]; 
title(tit); 
%%
plotsolu(s1,sol.y(:,1),12,2,2); title('q at t=0'); axis image; colormap hot; 
%% 
s1.plot.pstyle=3; mov=sol2mov(s1,sol,1);
%%
save('h2-17toFSC','sol','fn','alv');
