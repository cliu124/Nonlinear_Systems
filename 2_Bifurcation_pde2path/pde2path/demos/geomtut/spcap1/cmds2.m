%% cmc spherical cap, comparison of refinements; 
close all; keep pphome; 
global p2pglob; p2pglob.tsw=1; p2pglob.vi=[20, 40]; p2pglob.edc='k'; % plot-pars
%% tests of mesh-refinement, preparation: set rlong=1 and dsmax,dlammax 
p=loadp('cap1','pt10'); p.sw.rlong=1;  p.nc.dsmax=4; p.nc.dlammax=4; p0=p; 
%% alternate refine and cont, here ref.each 15th step; single steps for saving 
p=p0; p=setfn(p,'capr1'); sig=0.3; nsteps=13; 
for i=1:5; p=refineX(p,sig); p=cont(p,1); p=cont(p,1); p=cont(p,nsteps); end 
%% refine when error exceeds errbound, using refufu 
p=p0; p=setfn(p,'capr3'); p.fuha.ufu=@refufu; p.nc.errbound=0.04; p=cont(p,100); 
%% refine when max A exceeds p.maxA, using refufumaxA 
p=p0; p=setfn(p,'capr3'); p.fuha.ufu=@refufumaxA; p.fuha.e2rs=@e2rsmaxA; 
p.maxA=0.3; p=cont(p,100); 
%% error plots; error appended at end of cmcbra, component c=13 
lab=[25 27]; c=13; mclf(8); plotbra('capr1','pt71',8,c,'lab',lab,'fp',11); 
plotbra('capr3',8,c,'lab',[],'cl','r','fp',11); 
plotbra('capr4',8,c,'lab',[],'cl','m','fp',11); 
xlabel('V'); ylabel('||H-H(V)||_2/|H(V)|'); box on; 
axis([0,200,0.01,0.1]); 
%% h/r (mesh-quality), component c=npar+4=3+4=7
lab=[25 27]; c=7; mclf(8); plotbra('capr1','pt71',8,c,'lab',lab,'fp',11); 
plotbra('capr3',8,c,'lab',[],'cl','r','fp',11); 
plotbra('capr4',8,c,'lab',[],'cl','m','fp',11); 
xlabel('V'); ylabel('max(h/r)'); box on; 
%% Jac checks (timing and relative errors) 
for i=[40 42 72 102 133]; 
  p=loadp('capr2',['pt' mat2str(i,3)]); p.file.count=p.file.count-1; p.np, pplot(p,1); p.nc.del=0.0001; 
  p.branch(6+4,end), p.branch(6+5,end), p.branch(6+5,end)/p.np, p.branch(6+19,end)
 % [Gu,Gn]=jaccheck(p,0.5); [qu,qun]=qjaccheck(p,0.5); % slow!
  pause 
end
%% soln plots 
pplot('capr1','pt25'); pause; pplot('capr1','pt27'); pause, pplot('capr1','pt70');
%% alternate refine and cont, here ref.each 20th step; single steps for saving 
p=loadp('cap1','pt10','capr2'); p.sw.nobdref=0; p.sw.rlong=1; sig=0.3; nsteps=18; 
p.nc.dsmax=4; p.nc.dlammax=4;
for i=1:4; p=refineX(p,sig); p=cont(p,1); p=cont(p,1); p=cont(p,nsteps); end 