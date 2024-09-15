%% Brusselator auxcmds, illustrating Hopf-bif after mesh-refinement, movies and spectra
% 1D: meshrefine 1ds1/bpt2 and then cont 
p=loadp('1ds1','pt5','1dref'); p=meshada(p); p=cont(p,10);
%% Hopf bif on mesh-refined Turing branch 
ds=0.2; para=4; aux=[]; aux.tl=20; 
p=hoswibra('1dref','hpt2',ds,para,'1ds1h2ref',aux); 
p.nc.dsmax=1; p.sw.verb=2; p.hopf.jac=1; p=setbel(p,1,1e-4,5,@lss);
tic; p=cont(p,5); toc
%%  same in 2D 
p=loadp('2ds1','pt5','2ds1ref'); p=meshada(p,'sig',0.2);  
p.sol.restart=1; p.sol.ds=-0.005; p.nc.dsmax=0.005; p=cont(p);
%% Hopf bif on mesh-refined Turing branch 
aux=[]; aux.tl=10; para=4; 
p=hoswibra('2ds1ref','hpt1',0.1,para,'2ds1h1ref',aux); 
p.nc.tol=1e-6; p.file.smod=1; p=setbel(p,1,1e-4,5,@lss); p.hopf.flcheck=0; 
p.hopf.xi=1e-3;p.nc.dsmax=0.1; p.sw.verb=2; tic; p=cont(p,5); toc
%% a movie 
p=loadp('2ds1h1ref','pt5');  homov(p,1,1);
%% to illustrate problems with many small real eigenvalues, and need for initeig 
% modify bruinit to ly=lx! 
ndim=2; dir='hom2dbig'; p=[]; nx=100;
Du=0.01; Dv=0.1; Dw=1; c=1; d=1; a=0.95; b=2.75; lx=pi/1.4;  %lx=pi
par=[a b c d Du Dv Dw]; % a, b, c, d, Du, Dv, Dw 
p=bruinit(p,lx,nx,par,ndim); p=setfn(p,dir); 
%% compute many EVs close to 0 
r=resi(p,p.u); Gu=getGu(p,p.u,r); p.nc.neig=200; p.sw.verb=2; tic; ineg=spcalc(Gu,p); toc
%% compute om-guess 
p.wn.Msw=1; p=initwn(p,1,2); p=initeig(p,10); p.nc.neigv=[3, 3]; p0=p; 
figure(8); axis tight; set(gca,'fontsize',16); xlabel('i\omega'); ylabel('g(i\omega)'); 
%% compute evals near om-guesses 
p.nc.neigv=[3, 3];
[ineg,muv1]=spcalc(Gu,p,1); [ineg,muv2]=spcalc(Gu,p,2); pause
figure(6); clf; plot(real(muv1), imag(muv1),'b*'); hold on; 
plot(real(muv2), imag(muv2),'r*'); 