close all; keep pphome; % AC on sector 
%% continue trivial branch without rotational phase-condition (PC) 
p=[]; par=[-0.2]; phi=pi/6; R=5; nr=20; p=acinit(p,R,phi,nr,par); 
p=setfn(p,'trs'); p.nc.ilam=1; p.sw.bifcheck=2; p.nc.mu1=0.2; p=cont(p,50);
%%
for i=1:3; p.pdeo.grid.identifyBoundarySegment(i); 
    set(gca,'fontsize',14); axis image; pause; end 
%% bif branches 
ind=1:3; 
for i=ind
   is=mat2str(i); p=swibra('trs',['bpt' is],['b' is],0.1); p=cont(p,20); 
end
%% double the sector; first define dummy p on double sector and do 1 step for saving 
p=[]; par=[-0.2]; p=acinit(p,R,2*phi,nr,par); 
p=setfn(p,'dummy'); p.nc.ilam=1; p.sw.bifcheck=2; p.nc.mu1=0.5; p=cont(p,1);
%% reload 'dummy' (not strictly necessary, but useful for testing), load pt to double. 
p=loadp('dummy','pt1'); q=loadp('b1','pt10'); p=doublesec(p,q,phi); 
%% continue the doubled soln 
p=setfn(p,'b1d'); p=resetc(p); p=cont(p,20); 
%% BD plot 
fnr=3; figure(fnr); clf; cmp=4; 
plotbra('trs',fnr,cmp,'lsw',0);  plotbra('b1',fnr,cmp,'cl',p2pc('r1'),'lab',10);  
plotbra('b2',fnr,cmp,'cl',p2pc('r2'),'lab',10);  
plotbra('b3',fnr,cmp,'cl',p2pc('r3'),'lab',10); 
plotbra('b1d',fnr,cmp,'cl','m','lab',0); 
xlabel('\lambda'); ylabel('||u||_2'); 
%% soln plots 
plotsol('b1','pt10',1,1,2); nola; pause; plotsol('b2','pt10',1,1,2); nola; pause; 
plotsol('b3','pt10',1,1,2); nola; plotsol('b1d','pt0'); nola 
%% check mesh-adaption on original b1/pt10 and on doubled sector 
p=loadp('b1','pt10'); p.nc.maxt=2000; p=oomeshada(p,'ngen',2,'sig',0.2); pause 
%%
p=loadp('b1d','pt10'); p.nc.maxt=2000; p=oomeshada(p,'ngen',2,'sig',0.2); 