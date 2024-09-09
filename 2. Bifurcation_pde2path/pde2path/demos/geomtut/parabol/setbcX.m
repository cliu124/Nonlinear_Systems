function X=setbcX(p) % set BCs here, use u=0 on bdry in sG
X=p.X; id=p.idx; x=X(id,1); y=X(id,2); z=zfu(p,x,y); X(id,3)=z; 