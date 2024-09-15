function uf=nbcext(p,u)
nx=p.nx; ny=p.ny; uf=zeros((p.nx+2)*(p.ny+2),1);
uf(p.bui)=u(1:p.np); uf=reshape(uf,ny+2,nx+2); % pause 
uf(:,1)=uf(:,2); uf(:,nx+2)=uf(:,nx+1); 
uf(1,:)=uf(2,:); uf(ny+2,:)=uf(ny+1,:); 
uf=uf(:); 