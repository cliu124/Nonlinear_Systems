close all; format compact; keep pphome; 
%% init
ndim=4; dir='sqhom'; p=[]; lx=pi; par=[-0.2; 1; 0.1; -1; 1]; nx=10; 
p=cGLinit(p,lx,nx,par,ndim); p=setfn(p,dir); p.nc.mu1=1; p.nc.lammax=2; p.usrlam=[1 2]; 
p.sw.bifcheck=2; p=cont(p,10);
%% branch switching at 2nd BP which is double (hor.&vert stripes) 
para=4; ds=0.1; aux=[]; aux.tl=30; 
for sw=2; 
    switch sw; 
        case 1; aux.z=[1 0]; o='sq1'; % vertical stripes 
        case 2; aux.z=[1 -1]; o='sq2'; % vertex 
        case 3; aux.z=[1 -1i]; o='sq3'; % rot. osc. 
        case 4; aux.z=[0 1]; o='sq4'; % horiz. stripes 
    end 
    p=hoswibra('sqhom','hpt2',ds,para,o,aux); p.nc.lammax=2; p.usrlam=[1 2]; pause 
    p.hopf.jac=1; p.nc.dsmax=0.4; p.hopf.xi=1e-3; p.sw.verb=2; p.hopf.flcheck=0; 
    AMG=1;  % set AMG=1 if ilupack available
    if ~AMG; p=setbel(p,1,1e-3,10,@lss); % use BEL without ilupack 
    else  % use AMG with or without bel; AMG seems indifferent to borders! 
    p=setilup(p,1e-3,100); p.fuha.lss=@lss; p.fuha.blss=@lssAMG;
    end 
    tic; p=cont(p,30);  toc
end
%% 3rd HBP, simple 
p=hoswibra('sqhom','hpt3',ds,para,'sq5',aux); p.nc.lammax=2; p.usrlam=[1 2]; pause 
p.nc.dsmax=0.4; p.hopf.xi=1e-3; p.sw.verb=2; p.hopf.flcheck=0; 
if ~AMG; p=setbel(p,1,1e-3,10,@lss); % use BEL without ilupack 
else  p=setilup(p,1e-3,100); p.fuha.lss=@lss; p.fuha.blss=@lssAMG; end
tic; p=cont(p,30);  toc
%% BD, L^2
bpcmp=9; figure(3); clf; 
plotbra('sq1',3,bpcmp,'cl','b','lab',19); plotbra('sq2',3,bpcmp,'cl','k','lab',20); 
plotbra('sq3',3,bpcmp,'cl','r','lab',19); plotbra('sq5',3,bpcmp,'cl','m','lab',14); 
axis([-0.1 2 0 1.25]); 
xlabel('r'); ylabel('||u||_*'); 
%% soln plot
aux=[];  aux.cax='image'; aux.lay=[1 4]; aux.pind=1:5:16; 
hoplotf('sq1','pt19',1,1,aux); pause; hoplotf('sq2','pt20',1,1,aux); pause 
hoplotf('sq5','pt14',1,1,aux); pause 
aux.view=[100,60]; hoplotf('sq3','pt19',1,1,aux); 
%%
p=loadp('sq4','pt24'); p.plot.pstyle=2; p.hopf.ax='unif'; p.plot.cm='cool'; 
aux.lay=[1 4]; %aux.pind=1:5:16, 
hoplot(p,1,1,aux); % hoplot still to update !
%% BD-T
bpcmp=6; figure(3); clf; 
plotbra('sq1',3,bpcmp,'cl','b','lab',19); plotbra('sq2',3,bpcmp,'cl','k','lab',20); 
plotbra('sq3',3,bpcmp,'cl','r','lab',19); %plotbra('sq4',3,bpcmp,'cl','m','lab',24); 
%axis([1 2.6 0 1]); 
xlabel('r'); ylabel('T');