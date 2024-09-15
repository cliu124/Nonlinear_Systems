%% go branch backwards, cont-coarsening loop
p=loadp('capr1','pt56','capr1b'); p.sol.ds=-3.2; p=resetc(p);p.fuha.e2rs=@e2rsAi;
p.file.smod=5; p0=p;  sig=0.5; for i=1:8; p=coarsenX(p,sig); p=cont(p,5);   end 
%% coarsen via coarsufu
p=p0;p=setfn(p,'capr1d');p.fuha.ufu=@coarsufu; p.nc.delbound=8; p=cont(p,40); 
%% branch plot, (un)comment to select compo 
c=7; fn=8; ylab='\delta_{mesh}=max(h/r)'; %c=8; fn=8; ylab='a_{max}'; 
lab=[0 35]; mclf(fn); plotbra('capr1b','pt35',fn,c,'cl','k','lab',lab,'tyun','-*','lp',40); 
plotbra('capr1d','pt35',fn,c,'cl','m','lab',[],'lp',40,'tyun','-*'); 
xlabel('V'); ylabel(ylab); box on; grid on; 
%% soln plots 
p2pglob.tsw=1; p2pglob.vi=[20, 40]; p2pglob.edc='k'; p2pglob.cm='parula'; p2pglob.cut=0; 
pplot('capr1b','pt35'); pause; pplot('capr1d','pt35');
%% MCF, with initial large V; to handle meshing, alternate flow and coarsening 
% this may require trial and error to balance dt, flow-length nf, and 
% coarsening sigc.     First some graphics settings, then load and prepare:
p2pglob.cut=0; p2pglob.vi=[30,40]; p2pglob.cm='spring'; p2pglob.tsw=4; 
p=loadp('cap1r','pt15','mcf'); sigc=0.1;  dt=0.0005; nf=500; nplot=100; 
p.sw.nobdcoarsen=0; p.t=0; plotHK(p); figure(1); title('t=0');  % prepare MCF 
p.fuha.flowf=@mcff; t=0; ts=[]; p.fuha.e2rs=@e2rsAi; 
%% the MCF/coarsening loop; repeat this cell as desired 
for i=1:4; [p.X,t,ts]=geomflow(p,t,ts,dt,nf,nplot); p=coarsenX(p,sigc); end 
%% time series plot, A and V 
mclf(11); plot(ts(1,:),ts(2,:),ts(1,:), ts(3,:)); legend('A','V'); xlabel('t'); axis tight; 
set(gca,'fontsize',16);  pause; 
%% time series delta 
mclf(12); plot(ts(1,:),ts(4,:)); legend('\delta'); set(gca,'fontsize',14); grid on; 
%% other tests: ------go branch backwards, cont-degcoarsenX loop (doesn't quite work) 
%p=p0; p=setfn(p,'capr1c'); sig=0.2; for i=1:8; p=cont(p,5); p=degcoarsenX(p,sig,1,10); end 
%% a point for which coarsening in MCF doesn't quite work 
p=loadp('cap1r','pt22','mcf'); sigc=0.075; dt=0.00025; nf=600; nplot=200; 
p.fuha.flowf=@mcff; p.fuha.e2rs=@e2rsAi; 
%% sometimes degcoarsenX and/or moveX (instead of just coarsenX) are better 
sigc=0.2; cit=30; p.fuha.flowf=@mcff; t=0; ts=[];
for i=1:10; [p.X,t,ts]=geomflow(p,t,ts,dt,nf,nplot); p=degcoarsenX(p,sigc,1,cit); end 
%%
sigc=0.2; 
for i=1:4; 
    for j=1:2; [p.X,t,ts]=geomflow(p,t,ts,dt,nf,nplot); p=coarsenX(p,sigc); 
    end
    p=moveX(p,0.001,1); 
end 