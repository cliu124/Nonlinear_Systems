%% convergence and timing tests on 2db1 at r=2; 
para=4; ds=0.1; aux=[]; aux.tl=30; figure(2); clf; nsteps=20; AMG=1; 
for m=[20 30,40] 
    aux.tl=m; %p.file.smod=1; 
    p=hoswibra('hom2d','hpt1',ds,para,['m' mat2str(m)],aux); 
    p.hopf.jac=1; p.nc.dsmax=0.4; p.hopf.xi=1e-3; p.sw.verb=2; p.hopf.flcheck=1; 
    if ~AMG; p=setbel(p,1,1e-3,10,@lss); % use BEL without ilupack 
    else  % use AMG with or without bel; AMG seems indifferent to borders! 
     p=setilup(p,1e-3,100); p.fuha.lss=@lss; p.fuha.blss=@lssAMG;
    end 
    p.usrlam=2; pause; tic; p=cont(p,nsteps); toc
    p=cont(p,5); % additional steps to reach r=2 and beyond 
end 
%% BD, L^2
bpcmp=9; figure(3); clf; plotbra('m20','pt26',3,bpcmp,'lab',20,'cl','k'); 
plotbra('2m20',3,bpcmp,'lab',15,'cl','b'); 
axis([1 3 0 0.75]); xlabel('r'); ylabel('||u||_*'); 
%% soln plot
aux=[];  aux.pstyle=3; aux.xtics=[-2 2]; aux.ytics=[-1 1]; aux.lay=[1 4]; 
aux.pind=[1 5 9 14]; hoplotf('2m20','pt15',1,1,aux); 
%% Floq-timing
[muv1,muv2,ind]=floqap('m40','pt20'); 
%% differences in periods for different m 
p=loadp('m20','pt20'); getaux(p)', T20=p.hopf.T; 
p=loadp('m30','pt24'); getaux(p)', T30=p.hopf.T; 
p=loadp('m40','pt27'); getaux(p)', T40=p.hopf.T; 
T30-T20, T40-T30

