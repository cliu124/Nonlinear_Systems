close all; keep pphome; % AC on disk 
%% continue trivial branch without PC 
p=[]; par=[-0.2 0 0]; p=acinit(p,5,20,par); p=setfn(p,'tr'); 
p.nc.ilam=1; p.sw.bifcheck=2; p.nc.mu1=0.5; p=cont(p,30); 
%% initial cont, with small Eval
ind=1:1:4; 
for i=ind
   is=mat2str(i); p=swibra('tr',['bpt' is],['b' is],0.1); p.sw.verb=2; pause; p=cont(p,2); 
end
%% cont 1st and 4th branch further, and other branches with rot.PC
p=loadp('b1','pt2'); p=cont(p,10); p=loadp('b4','pt2'); p=cont(p,10); 
for i=ind(2:end-1)
  is=mat2str(i); try; p=loadp(['b' is],'pt2'); catch; p=loadp(['b' is],'pt3'); end
  p.file.smod=10; p.nc.ilam=[1;2];  p.nc.nq=1; 
  p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; % analytical jac for aux. eqn.
  p.sw.bprint=[1 2]; clf(2); p.nc.dsmax=0.11;  p=cont(p,9); 
end  
%% BD plot (lam)
figure(3); clf; plotbra('tr','bpt4','lsw',0);  plotbra('b1','lsw',0,'cl',p2pc('b1'));  
plotbra('b2','lsw',0,'cl',p2pc('r1'));  plotbra('b3','lsw',0,'cl',p2pc('r2'));  
plotbra('b4','lsw',0,'cl',p2pc('b2'),'ms',0);  xlabel('\lambda'); ylabel('||u||_2'); 
%%
plotsol('b2','pt10',1,1,2); nolab; pause; plotsol('b3','pt11',1,1,2); nolab
%%
plotsol('b4','pt13',1,1,1); view(-20,60); zticks([0 0.5]); nolab 
%% plot lap(u)
p=loadp('b4','pt13','t1'); Ms=p.mat.M; Ks=p.mat.K; 
ulap=-Ms\(Ks*p.u(1:p.nu)); ulap2=-Ms\(Ks*ulap(1:p.nu)); 
plotsolu(p,ulap,6,1,1); 