%% AC on (small) sphere: eigenfunctions spherical harmonics, dim(ker)=2l+1, 
% with eigenfunctions pairwise related by spatial shifts (sin(mx) and cos(mx)). 
% hence, use aux.ali to reduce dim(ker) to l+1 in qswibra(even l) and cswibra (odd l) 
close all; keep pphome; 
%% init and trivial branch 
p=[]; par=[3 -0.1 1 0]; % parameters [R lambda gamma s]  (s=speed for x-PC) 
lx=pi; del=1e-3; ly=pi/2-del; nx=31; ny=17; ref=1;% nx=11; ny=5; ref=1;
p=acinit(p,lx,ly,nx,ny,par,ref); p=setfn(p,'tr3'); p=cont(p,20);
%% BP1, simple, spat.homogen.branch(es) 
p=swibra('tr3','bpt1','1a',0.1); p=cont(p,20); 
p=swibra('tr3','bpt1','1b',-0.1); p=cont(p,20); 
%% BP2, l=1, pitch, dim(ker)=2l+1=3, use aux.ali to select kernel vectors
aux=[]; aux.besw=0; aux.m=3;  aux.isotol=1e-3;
aux.ali=[1 3]; aux.besw=1; % in first run, comment out this line to see EFus
p0=cswibra('tr3','bpt2',aux); p0.sw.bifcheck=0; p0.nc.tol=1e-6; 
%% one spot
p=seltau(p0,2,'2-1',3); p.file.smod=5; p=conpc(p,0.05,2,20); 
%% BP3, l=2, dim(ker)=2l+1=5, trans, 2 branches (modulo symmetry) 
aux=[]; aux.isotol=1e-2; aux.besw=0; aux.m=5;  
aux.ali=[1 3 5]; aux.besw=1; % selecting kernel vectors 
p0=qswibra('tr3','bpt3',aux); 
p0.sw.bifcheck=0; p0.nc.tol=1e-6; 
%% stripe like 
p=seltau(p0,1,'3-1a',2); p=conpc(p,0.1,3,20); spplot(p); 
%p=seltau(p0,1,'3-1b',2); p.sol.ds=-0.04; p=conpc(p,0.1,3,20); spplot(p); 
%% horizontal stripe, hence no PC 
p=gentau(p0,[0 0 1],'3-2a'); p=cont(p,20); spplot(p);  % other direction falls onto 0
%% BP4, l=3, pitchforks 
aux=[]; aux.besw=0; aux.m=8; aux.isotol=0.05; % decrease isotol due to close by solns 
aux.ali=[1 3 5 7]; aux.besw=1;
p0=cswibra('tr3','bpt4',aux); p0.sw.bifcheck=0; p0.nc.tol=1e-6; 
%% 4 spots, tetra-sym 
p=seltau(p0,1,'4-1',3); p=conpc(p,0.02,3,15); 
%% spot+stripe
p=gentau(p0,[0 0 0 1],'4-3'); p=cont(p,20); spplot(p); 
%% plot BD 
f=3; c=0; figure(f); clf;
plotbra('tr3',f,c,'cl','k','lsw',0); 
plotbra('1a',f,c,'cl','m','lsw',0); plotbra('1b',f,c,'cl','m','lsw',0) ;
plotbra('2-1',f,c,'cl',p2pc('b1'),'lab',[3 14]); 
plotbra('3-1a',f,c,'cl',p2pc('r1'), 'lab',11); %plotbra('3-1b',f,c,'cl',p2pc('r1'),'lab',13) ;
plotbra('3-2a',f,c,'cl',p2pc('r3'),'lab',13); 
plotbra('4-1',f,c,'cl',p2pc('o1'),'lab',14);
%plotbra('4-2',f,c,'cl',p2pc('o2'),'lab',10); 
plotbra('4-3',f,c,'cl',p2pc('o3'),'lab',15);
axis([-0.5 3.1 0 5]); xlabel('\lambda'); ylabel('||u||_2'); 
%% soln plots 
spplotf('2-1','pt2'); pause; spplotf('2-1','pt14'); pause; spplotf('3-1a','pt11'); pause; 
spplotf('3-1b','pt13'); pause; spplotf('3-2a','pt13'); pause; 
%%
spplotf('4-1','pt14'); pause; spplotf('4-2','pt10'); pause; spplotf('4-3','pt15'); 
%%
[Ga,Gn]=jaccheck(p); Gd=abs(Ga-Gn); e1=max(max(Gd)); figure(10); clf; 
spy(Gd>e1/10); 
%% c8: continue in some other param, here the radius R
p=swiparf('3-1a','pt10','b3R',[1 4]); p.sol.ds=0.1; 
p.nc.lammin=-1; p.nc.lammax=20; clf(2); p=cont(p,20); 
%%
spplotf('b3R','pt1'); pause; spplotf('b3R','pt10'); pause; spplotf('b3R','pt20'); 