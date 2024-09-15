function [Ja,Jn]=hpjaccheck(p)
% HPJACCHECK: compare numerical and assembled pa_u(G_u phi) (for HP cont) 
u=p.u; p=hpreduce(p); nu=p.nu; % set regular case sizes and turn off spcont in this block
tic; Ja=p.fuha.spjac(p,u); t1=toc; % fast way 
figure(p.plot.spfig); spy(Ja); title('user spjac'); 
fprintf('time for assembling=%g\n',t1);
u1=[u(1:p.nu); u(3*p.nu+1:3*p.nu+p.naux)]; 
phir=u(p.nu+1:2*p.nu); phii=u(2*p.nu+1:3*p.nu);
nu1=length(u1); r1=pderesi(p,u1); 
Gv=getGupde(p,u1,r1); % get pde part linearization
Gvphr=Gv*phir; Gvphi=Gv*phii; 
tic; 
if p.sw.spjac==1; spjac=0; else spjac=p.sw.spjac; end 
if spjac==0;  Gvvph=getGuphih(p,u1,phir,phii); % numjac 
else 
Gvvph=sparse(2*nu,nu); del=p.nc.del; 
for j=1:nu  
  up=u1+del*ej(j,nu1); rp=pderesi(p,up); 
  Gvp=getGupde(p,up,rp); % perturbed pde-part lin.
  rr=(Gvp*phir-Gvphr)/del; ri=(Gvp*phii-Gvphi)/del;   
  Gvvph=Gvvph+sparse(1:2*nu,j,[rr;ri],2*nu,nu); 
end  
end
t2=toc; Jn=Gvvph; 
figure(p.plot.ifig); clf; spy(Jn); title('numerical hpjac'); 
fprintf('time for FD appr=%g\n',t2);
m1=full(max(max(abs(Jn)))); m2=full(max(max(abs(Jn-Ja)))); 
m3=full(max(sum(abs(Jn)))); m4=full(max(sum(abs(Jn-Ja))));
fprintf('max(Jn)=%g, max(Jn-Ja)=%g, infinity-norm(Jn)=%g, relerr=%g\n',...
    m1,m2,m3,m4/m3);
