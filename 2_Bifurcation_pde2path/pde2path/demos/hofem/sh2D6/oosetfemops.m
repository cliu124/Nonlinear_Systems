function p=oosetfemops(p) % for SH as 2nd order system, hence singular p.mat.M  
if p.hofem.sw;  % 6-node, simple syntax (c=a=1) [K,M]=assem6(p);
 try K=load(p.hofem.Kfn,'KM'); K=K.KM.K; M=K.KM.M; % try loading K,M from disk 
% if that fails, assemble K,M and save 
 catch;tic; [K,M]=assem6(p); toc, KM.K=K; KM.M=M; save(p.hofem.Kfn,'KM'); end
else [K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % standard 3-node
end 
p.mat.M=[[M 0*M];[0*M 0*M]];  % system mass matrix (here singular)
p.mat.Ks=K; p.mat.Ms=M;   % save SCALAR Laplacian, system K composed in sG