%% meander from iguess 
lx=pi/2; ly=3*pi/4; p.pp=3; p.del=0; sw=3; wi=0.25; a=1.25; 
nx=40; ny=40; p.fuha.wfu=@wfu2; % potential and derivatives 
p.eta1=1; p.eta2=3; p.eps=0.35; p.mi=-0.6; 
p=fchinits(p,lx,ly,nx,ny,wi,a,sw); p=setfn(p,'ms'); 
%% try to converge! pick suitable mesh by hand!
p=p0; 
p.fsol.fsol=1; % only works if first 3 refinement steps with fsolve as solver
p.nc.ngen=1; p.fsol.imax=5; p=meshref(p,'maxt',30000);
p.fsol.fsol=0; p.nc.ngen=1;
p=meshref(p,'maxt',30000);p=setbmesh(p);p.fuha.savefu(p);
%% continue ms branch 
p.nc.lammax=2; p=pmcont(p,20);
%% bifurcate 
p=swibra('ms','bpt1','ms1',-0.02); p.nc.dsmax=0.04; p=pmcont(p,20);
p=swibra('ms','bpt2','ms2',-0.02); p.nc.dsmax=0.04; p=pmcont(p,20);
%% plot BD 
figure(3);clf;cmp=3;
plotbra('ms',3,cmp,'cl','k','ms',0,'lwun',3,'lwst',3);
plotbra('ms1',3,cmp,'cl','b','ms',0,'lwun',3,'lwst',3);
plotbra('ms2',3,cmp,'cl','r','ms',0,'lwun',3,'lwst',3);
xlabel('\eta_1'); ylabel('\gamma'); title(['mass=-0.6']);
%% sol plots 
plotsol('ms','bpt1',1,1,2); xlabel(''); ylabel(''); pause 
plotsol('ms1','pt10',1,1,2); xlabel(''); ylabel(''); pause
plotsol('ms2','pt10',1,1,2); xlabel(''); ylabel(''); 
%%
del=0.05; newtontest('ms1','pt10',del); 