%% cgl, 3D, DBC, expensive! 
close all;  keep pphome; 
%% init and find HBPs 
ndim=3; dir='hom3d'; p=[]; lx=pi; nx=25; 
par=[5; 1; 0.1; -1; 1]; %  r  nu  mu   c3  c5
p=cGLinit(p,lx,nx,par,ndim); p=setfn(p,dir); p.sw.bifcheck=2; 
p.nc.mu1=0.5; p.plot.pstyle=1; p=cont(p,20);
%% hoswibra and cont of Hopf branches; this very much needs ILUPACK
para=4; ds=0.2; aux=[]; aux.tl=20; nsteps=15; figure(2); clf;
for bnr=1:2
switch bnr
 case 1; p=hoswibra('hom3d','hpt1',ds,para,'3db1',aux); 
 case 2; p=hoswibra('hom3d','hpt2',ds,para,'3db2',aux);  
end 
p.hopf.jac=1; p.nc.dsmax=0.85; p.hopf.xi=5e-3; 
p.hopf.flcheck=0; % switch off floquet (somewhat slow) 
p.file.smod=5; p.sw.verb=2; % reset this for, e.g., more output 
AMG=1; p.sw.verb=2; % set AMG=1 if ilupack available (if no ilupack, bel might give memory problems!) 
if ~AMG; p=setbel(p,1,1e-4,20,@lss); % use BEL without ilupack 
else  % use AMG with or without bel; AMG seems indifferent to borders! 
    %p=setbel(p,0,beltol,belimax,@lssAMG); p=setilup(p,droptol,AMGmaxit);
    p=setilup(p,1e-3,100); p.fuha.lss=@lss; p.fuha.blss=@lssAMG;
end 
t1=tic; p=cont(p,nsteps); toc(t1) 
end 
%% plot BD, L^2
bpcmp=9; figure(3); clf; plotbra('3db1','pt15',3,bpcmp); 
plotbra('3db2',3,bpcmp,'cl','b','lab',15); 
%axis([5.25 6.3 0 0.72]); set(gca,'XTick',[5.5 6]);
xlabel('r'); ylabel('||u||_*');
%% plot BD, T
bpcmp=6; pstyle=3; figure(3); clf; plotbra('3db1',3,bpcmp); 
plotbra('3db2',3,bpcmp,'lab',10, 'cl','b'); 
axis([5.5 7 6.4 7.2]); xlabel('r'); ylabel('T');
%% soln plot with custom layout (slice-plot at t=0, T/2) 
aux=[]; aux.lay=[1 2]; aux.pind=[1 10]; aux.pstyle=2; % slice-plot
aux.ztics=[-0.5 0.5]; hoplotf('3db2','pt10',1,1,aux); 
%% isoplot
aux.pstyle=0; hoplotf('3db2','pt10',1,1,aux); 
%% movie, in matlab play with movie(mov) 
p=loadp('3db2','pt15'); mov=homov3d(p,1,1,aux); 
mymov2avi(mov,'mcGL3d');
%% Floquet a posteriori, this is slow!
aux.nfloq=50; [muv1,~,~,~,h]=floqap('3db1','pt10',aux); 