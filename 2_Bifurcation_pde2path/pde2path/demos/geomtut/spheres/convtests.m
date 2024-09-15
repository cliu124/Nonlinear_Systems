%% experiments re convergence of discrete H and K 
p2pglob.cm='spring'; p2pglob.Hf=10; p2pglob.cut=7; p2pglob.showN=0; 
p2pglob.cb=1; p2pglob.vi=[60 60]; 
par=[0;0;0;0;0;0]; mclf(10); mclf(11); r0=1; sw=2; 
p=sphereinit(par,r0,sw); p.up=p.u; fn=['t' mat2str(sw)]; p=setfn(p,fn); plotHK(p); 
np=p.np, pause, N=getN(p,p.X); M=massmatrix(p.X,p.tri,'full'); % check full M 
LB=cotmatrix(p.X,p.tri); H2=0.5*dot(LB*p.X,N,2); H2=M\H2; % cot-H 
[k,H1,K1]=discrete_curvatures(p.X,p.tri); K2=M\K1; 
p.up(1:p.np)=H2; pplot(p,10); title(['H_{full}, np=' mat2str(np)]); colormap autumn; 
p.up(1:p.np)=K2; pplot(p,11); title(['Kf, np=' mat2str(np)]); colormap winter; 