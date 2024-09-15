%% Bratu model
close all; keep pphome; 
%% Homogeneous solution branch
p=[]; p=bratu_init(p,0.5,0.5,20); p=pmcont(p);  
%% 2 branch switchings and successive continuation
p=swibra('p','bpt1','q'); p.nc.lammin=0.09; p.sw.foldcheck=1; p=cont(p,20); 
p=swibra('q','bpt1','r'); p.nc.lammin=0.05;p=cont(p,20); 
%% Postprocessing, first plot BifDiagram 
figure(3);clf(3); plotbra('p',3,0,'cl','k'); 
plotbra('q',3,0,'lab',5,'cl','b'); plotbra('r',3,0,'lab',5,'cl','r');
axis([0.05 0.37 0.5 3.5]);xlabel('\lambda');ylabel('||u||_2');
%% Plot solutions 
plotsolf('q','pt5',7,1,1); % mesh plot 
figure(7); view(-30,50); xlabel('x');ylabel('y');
plotsolf('r','pt5',8,1,3); % rendered 3D plot 
figure(8); view(-30,50); xlabel('x');ylabel('y');
%% fold point continuation 
pb=spcontini('q','fpt1',2,'pb');    % init branch point continuation in par 3
pb.plot.bpcmp=pb.nc.ilam(2); clf(2); % use this new parameter for plotting
pb.nc.tol=1e-5;                  % increase tolerance as typically required for fold cont.
pb.sol.ds=1e-3;                   % new stepsize in new primary parameter
pb.sw.spjac=1; pb.fuha.spjac=@bratu_spjac; 
pb.nc.nsteps=20; pb.nc.amod=0; pb.sw.bifcheck=0; pb=cont(pb);
%% Exit fold point continuation and get tangent
clf(2); p2=spcontexit('pb','pt20','pt2'); p2.sol.ds=0.01; p2.tau=getinitau(p2);
%% Switch to nontrivial branch and continue
q2=swibra(p2,'q2',-0.05);
q2.nc.ilam=1; clf(2);
q2.plot.bpcmp=0; q2.nc.dsmin=1e-5;
q2.nc.amod=0; q2=cont(q2);
%% continuation with mesh-adaption depending on error bound 
q=swibra('p','bpt1','q'); q.nc.lammin=0.09;
q.sw.errcheck=2; q.nc.errbound=0.01; q=cont(q,20); 