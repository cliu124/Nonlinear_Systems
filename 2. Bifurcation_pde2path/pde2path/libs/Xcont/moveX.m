function p=moveX(p,deltat,it) 
% moveX:  remesh X and move points; taken from distmesh
% 
% p=moveX(p,deltat,it); 
X=p.X; t=p.tri; fh=@huniform; pfix=p.idx; h0=1; % adapting to distmesh-names 
dptol=1e-4; ttol=.1; Fscale=1.2; 
[t2t,t2n]=mkt2t(t); %  element connectivities
t2t=int32(t2t-1)'; t2n=int8(t2n-1)';
N=size(X,1);                                         % Number of points N
Xold=inf;          count=0;                          % For first iteration
while count<it
  X0=X; count=count+1; 
  % 3. Retriangulation
  if max(sqrt(sum((X-Xold).^2,2))/h0)>ttol           % Any large movement?
    Xold=X;                                          % Save current positions
    [t,t2t,t2n]=trisurfupd(int32(t-1)',t2t,t2n,X');  % Update triangles (mex, the crucial step!) 
    t=double(t+1)';
    % 4. Describe each bar by a unique pair of nodes
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
    bars=unique(sort(bars,2),'rows');                % Bars as node pairs
    % 5. Graphical output of the current mesh
    mclf(20),patch('faces',t,'vertices',X,'facecol',[.8,.9,1],'edgecol','k');
    axis equal; view(3); drawnow, %pause  
  end
  % 6. Move mesh points based on bar lengths L and forces F
  barvec=X(bars(:,1),:)-X(bars(:,2),:);              % List of bar vectors
  L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
  hbars=feval(fh,(X(bars(:,1),:)+X(bars(:,2),:))/2); %,varargin{:});
  L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
  F=max(L0-L,0);                                     % Bar forces (scalars)
  Fvec=F./L*[1,1,1].*barvec;                         % Bar forces (x,y,z components)
  Ftot=full(sparse(bars(:,[1,1,1,2,2,2]),ones(size(F))*[1,2,3,1,2,3],[Fvec,-Fvec],N,3));
  Ftot(pfix,:)=0; 
  X=X+deltat*Ftot;                                   % Update node positions
  % 8. Termination criterion: All nodes move less than dptol (scaled)
  if max(sqrt(sum((X-X0).^2,2))/h0)<dptol, break; end
end
p.X=X; p.tri=t; p.np=size(X,1); p.nt=size(t,1); 
end
