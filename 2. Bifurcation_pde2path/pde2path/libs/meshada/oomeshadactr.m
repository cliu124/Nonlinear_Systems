function p=oomeshadactr(p,varargin)
% OOMESHADACTR: trullekrul version
% important (hidden) parameters: p.trcop.npb=max np goal for coarsening step 
% trop=(essentially) trullekrul options (+npb and sw=control for adaptm) 
%
% called by oomeshadac if p.sw.trul>0 
if isfield(p,'trcop'); 
trops=p.trop; p.trop=p.trcop; 
np=p.np; npo=np; %p.trop.sw=4; p.trop.Llow=1000; options.fastRM=1; 
count=0; try; crmax=p.trcop.crmax; catch crmax=3; end 
while np>p.trop.npb && count<crmax; 
    p=oomeshada(p,'ngen',1); count=count+1; 
    np=p.np;
end
fprintf('old mesh coarsened from np=%i to np=%i\n',npo,np); 
p.trop=trops; % reset trop for refinement 
end 
p=oomeshada(p,varargin{:}); % refine the coarse mesh