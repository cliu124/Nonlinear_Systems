function p=oosetfemops(p) % for SH as 2nd order system, hence singular p.mat.M  
if p.hofem.sw % 10-node tetras, assembly slow, hence try loading from disk 
  try KM=load(p.hofem.Kfn,'KM'); K=KM.KM.K; M=KM.KM.M; % try to load K,M from disk 
  catch;tic; [K,M]=assem10(p); toc % assemble 
  KM.K=K; KM.M=M; save(p.hofem.Kfn,'KM');   % and save 
  end 
else [K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % OOPDE assembly 
end
p.mat.M=[[M 0*M];[0*M 0*M]];  % system mass matrix (here singular)
p.mat.Ks=K; p.mat.Ms=M;   % SCALAR Laplacian, system K composed in sG