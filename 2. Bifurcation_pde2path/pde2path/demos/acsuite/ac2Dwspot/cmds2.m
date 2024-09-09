%% continue soln at xi=0 in lambda, swiparf, then reset some parameters
p=swiparf('wsada','pt16','lamc',2); p.nc.dsmax=0.01; p.sol.ds=0.01; 
p.sw.foldcheck=0; p.nc.mu2=0.01; p.nc.foldtol=0.2; % reset foldtol (poor localization) 
p.trcop.npb=2000; p.trop.etafac=5e-6; p=cont(p,100); % allow finer meshes, and go 
%% 1st Bif. branch 
p=swibra('lamc','bpt1','lc1',-0.01); pause; p.nc.dsmax=0.025; p.sw.bifcheck=0; 
p.usrlam=[]; p.nc.amod=20; p.trcop.crmax=1; p.nc.neig=20; p=cont(p,50); 
%% plot BD
figure(3); clf;
plotbra('lamc','pt40',3,0,'cl','k','lab',[20 40 60],'fms',0,'lp',63,'fp',4);
plotbra('lc1','pt50',3,0,'cl','r','lab',[20 30 40 50],'lp',55);
xlabel('\lambda'); ylabel('||u||_2'); 
%% soln plots, lamc
iv=10:10:60; va=[50 60]; 
for i=iv 
    p=loadp('lamc',['pt' mat2str(i)]); 
    plotsol(p); myticks(['pt' mat2str(i) ', n_p=' mat2str(p.np)]); 
    view(va); set(gca,'fontsize',16); 
    p.np, pause; 
end 
%% soln plots, lc1
iv=20:10:50; va=[50 60]; 
for i=iv 
    p=loadp('lc1',['pt' mat2str(i)]); 
    plotsol(p); myticks(['pt' mat2str(i) ', n_p=' mat2str(p.np)]); 
    view(va); set(gca,'fontsize',16); 
    p.np, pause; 
end 
