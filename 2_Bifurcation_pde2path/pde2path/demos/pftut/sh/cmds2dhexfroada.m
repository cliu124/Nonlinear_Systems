%% refinement of hex2zero front by trullekrul and subsequent cont with amod=5
p=loadp('2D/H8f','pt60','2D/hexfroref'); p.nc.tol=1e-8; p=resetc(p); p.np 
op=troptions2D(); % load default trullerup-options, then overload some 
op.innerit=1; op.etafu=@etafua2D; op.verbose=2; op.ppar=2; op.Lup=2; 
op.Llow=0.075; % small LLow important to avoid too coarse meshes 
p.trop=op; ops=op; % put options in p, then modify to coarsening options 
op.npb=3000; op.sw=5; op.innerit=3; op.crmax=2; p.trcop=op; p.sw.trul=1; 
p.sw.ips=2; % interpolation by 'nearest', otherwise bad behaviour at boundary  
% call pure coarsening, then reset options 
p.trop=p.trcop; p=oomeshada(p,'ngen',3); p.trop=ops; stansavefu(p); 
%% meshada each 5th step, switch off pure coarsen(not needed/effective) 
p=loadp('2D/hexfroref','pt0','2D/hexfroada'); huclean(p); p.nc.tol=1e-6; p.trcop.crmax=0; 
p.nc.amod=5; p.nc.dsmax=0.01; p=resetc(p); stansavefu(p); 
p.sw.foldcheck=0; p.trop.innerit=2; p.nc.ngen=2; p=cont(p,101); 
%% BD 
figure(3); clf
plotbra('2D/hexfroada','pt100',3,3,'cl','m','lab',[0,10,80]); 
plotbra('2D/H8f',3,3,'fp',60,'cl','r','fms',0); 
xlabel('\lambda'); ylabel('||u_1||_*'); 
axis([-0.2 -0.175 0.19 0.35]); 
%% Sol plots
iv=[0 10 80]; 
for i=iv 
  p=loadp('2D/hexfroada',['pt' mat2str(i)]); p.np, 
  plotsol(p,11,1,2); xlabel(''); ylabel(''); 
  %plotsol(p,12,1,0); axis([-lx 0 -ly 0]);  xlabel(''); ylabel(''); 
  plotsol(p,12,1,1); view(5,20); zticks([]); axis([0 lx -ly 0]);  
  xlabel(''); ylabel(''); title(['hexf-ada, pt' mat2str(i) ', n_p=' mat2str(p.np)]); 
  pause
end 
