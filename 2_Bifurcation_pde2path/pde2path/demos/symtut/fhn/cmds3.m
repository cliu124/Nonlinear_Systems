%% cell 1 - use p3 to break symmetry to obtain direction reversing front
p=swiparf('stand','pt7','asym',[3;2]); p.nc.lammin=-0.4; p.nc.lammax=1; clf(2); 
p.usrlam=0.05; p.sol.ds=0.1; p.sw.foldcheck=1; p.plot.bpcmp=2; p=cont(p,10);
%% cell 2 - continue the fold to get the underlying cusp
p=spcontini('asym','fpt1',4,'cusp'); p.plot.bpcmp=p.nc.ilam(2); clf(2); 
p.sw.spjac=1; p.fuha.spjac=@spjac; p.sw.spqjac=1; p.fuha.spqjac=@spqjac; 
p.sol.ds=-0.01; p.nc.lammax=1.5; p=cont(p); 
%% cell 3 - perturb by eigenfunction and time-integrate 
p=loadp('asym','pt2','sim'); [muv,V]=specGu(p); n=p.nu; 
p.u(1:n)=p.u(1:n)+1e-3*real(V(1:n,1)); vel1=[]; t1=0; nt=5000; dt=0.5; pmod=500; vmod=50; 
[q,t1,vel1]=tintfreeze(p,t1,dt,nt,pmod,vel1,vmod);
p.u(1:n)=p.u(1:n)-2e-3*real(V(1:n,1)); vel2=[]; t1=0; % other direction 
[p,t1,vel2]=tintfreeze(p,t1,dt,nt,pmod,vel2,vmod);
p=setaux(p,2,vel2(2,end)); p.fuha.savefu(p); 
%% cell 4 - velocity and position plots
figure(7); clf; fs='fontsize'; 
plot(vel1(1,2:end),vel1(2,2:end),vel2(1,2:end),vel2(2,2:end),'-r','LineWidth',3); 
xlabel('time',fs,16); ylabel('velocity',fs,16); set(gca,fs,16); axis tight
tl=size(vel1,2); pos1=zeros(1,tl); pos2=pos1; 
for(i=1:tl-1) 
    pos1(i+1)=pos1(i)+vel1(2,i)*(vel1(1,i+1)-vel1(1,i));
    pos2(i+1)=pos2(i)+vel2(2,i)*(vel2(1,i+1)-vel2(1,i));
end;
figure(8); clf; plot(vel1(1,:),pos1,vel2(1,:),pos2,'-r','LineWidth',3)
xlabel('time',fs,16); ylabel('position',fs,16); set(gca,fs,16); axis tight
%% cell 5 - continue final state as TW (get branch from cell 1 again!)
p=swiparf('sim','pt3','travel2',[3;2]); p.usrlam=p.u(p.nu+3); p.nc.lammin=-4; 
p.nc.lammax=0.2; p.sol.ds=0.1; p.nc.dsmax=0.1; clf(2); p=cont(p,30);
figure(3); clf; plotbra(p,3,2,'lsw',0);