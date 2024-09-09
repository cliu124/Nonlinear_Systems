%% cell 1 - perturb by eigenfunction and simulate time dynamics
p=loadp('stand','pt7','sim'); p.nc.nq=0; [muv,V]=specGu(p); muv(1:4)
p.u(1:p.nu)=p.u(1:p.nu)-1e-3*real(V(1:p.nu,2)); 
vel1=[]; t1=0; nt=5000; dt=0.5; pmod=500; vmod=50; 
[q,t1,vel1]=tintfreeze(p,t1,dt,nt,pmod,vel1,vmod);
t1=0; vel2=[]; p.u(1:p.nu)=p.u(1:p.nu)+2e-3*real(V(1:p.nu,2)); % other direction: 
[q,t1,vel2]=tintfreeze(p,t1,dt,nt,pmod,vel2,vmod);
%% cell 2 - velocity plot
figure(7); clf; fs='fontsize'; 
plot(vel1(1,:),vel1(2,:),vel2(1,:),vel2(2,:),'-r','LineWidth',3); 
xlabel('time',fs,16); ylabel('velocity',fs,16); set(gca,fs,16); axis tight
%% cell 3 - position plot
tl=size(vel1,2); pos1=zeros(1,tl); pos2=pos1; 
for(i=1:tl-1) 
    pos1(i+1)=pos1(i)+vel1(2,i)*(vel1(1,i+1)-vel1(1,i));
    pos2(i+1)=pos2(i)+vel2(2,i)*(vel2(1,i+1)-vel2(1,i));
end;
figure(8); clf; plot(vel1(1,:),pos1,vel2(1,:),pos2,'-r','LineWidth',3)
xlabel('time',fs,16); ylabel('position',fs,16); set(gca,fs ,16);axis tight