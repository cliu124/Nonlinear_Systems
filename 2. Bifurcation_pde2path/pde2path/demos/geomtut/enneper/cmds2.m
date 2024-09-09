%% MCF; to handle meshing, alternate flow, degcoarsening  and moveX 
% this requires some trial and error, i.e., flow-length (ns) vs coarsening sig
p2pglob.cut=0; p2pglob.cm='spring'; p2pglob.cb=1; p2pglob.showbd=2; p2pglob.vi=[-10 25];
p=loadp('e1','pt12','dummy'); p.t=0; par=getaux(p); par(4) % check alpha
pXf=0.05*(rand(p.np,1)-0.4); pXf(p.idx)=0; N=getN(p,p.X); p.X=p.X+pXf.*N; 
p.fuha.flowf=@mcff; t=0; ts=[]; dt=0.002; ns=500; nplot=100; % prepare MCF 
sigr=0.025; p.sw.rlong=1; p.u(1:p.nu)=0; p2pglob.tsw=4; p.fuha.e2rs=@e2rsA; 
sigc=0.05; cit=10; keepbd=1; mdt=0.0001; mit=50; % parameters for coarsen and move 
clear mov; mc=0; % for movie 
%% flow, with refine, coarsen and moveX, and getframe for movie 
for m=1:4   
for i=1:4; [p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot);
    if p.np<1100; p.fuha.e2rs=@e2rsA; p.sw.nobdref=0; 
     tho=p.th; idold=p.idx; p=refineX(p,sigr); p=thinterpol(p,idold,tho); 
     Xbc=bcX(p,p.u); p.X(p.idx,:)=Xbc; 
    end       
    p.t=t; pplot(p); mov(mc+1)=getframe(1); mc=mc+1; 
    if p.np>700; p=degcoarsenX(p,sigc,cit,keepbd); end     
end
p=moveX(p,mdt,mit);
end 
p1=p; ts1=ts; 
%% time series plot 
mclf(11); plot(ts(1,:),ts(2,:),ts(1,:), ts(3,:),'linewidth',2); 
legend('A','V'); axis tight; set(gca,'fontsize',14); xlabel('t'); 
%%
mymov2avi(mov,'m1'); 