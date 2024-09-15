function p=oosetfemops(p) % Schnak on cone 
par=p.u(p.nu+1:end); a=par(4); sw6=0; try sw6=p.sw6; catch; end 
if sw6; 
 try KM=load(p.Kfn,'KM'); K=KM.KM.K; M=KM.KM.M; % try to load K,M from disk 
  catch;tic; [K,M]=LBcone6(p,a); toc % assemble 
  KM.K=K; KM.M=M; save(p.Kfn,'KM');   % and save 
 end 
else [K,M]=LBcone(p,a); end 
M=[M 0*M; 0*M M]; p.mat.M=M; p.mat.K=K;
