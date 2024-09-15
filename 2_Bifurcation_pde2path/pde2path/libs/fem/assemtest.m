%
% assemtest: SCRIPT for testing assema, dummy, convenience, freely modify!
np=p.np; nt=p.pdeo.grid.nElements; fem=p.pdeo.fem; 
gr=p.pdeo.grid; po=getpte(p); x=po(1,:); y=po(2,:);
%% 2D
cc1=1;  [K1,M1]=fem.assema(gr,cc1,1,0); 
cc2=[1 0; 0 1]; [K2,M1]=fem.assema(gr,cc2,1,0); 
max(max(abs(K1))), max(max(abs(K2)))
cc3=[0 1; 1 0]; [K3,M1]=fem.assema(gr,cc3,1,0); 
ov1=ones(1,np); ov2=ones(1,nt); 
cc4=[1*ov1 0*ov1; 0*ov1 1*ov1]; [K4,M1]=fem.assema(gr,cc4,1,0); 
K1-K2, K3-K4, max(max(abs(K3))), max(max(abs(K4)))
cc5=[0*ov2 ov2; ov2 0*ov2]; [K5,M1]=fem.assema(gr,cc5,1,0); 
K4-K5
%% 3D 
ov1=ones(1,np); ov2=ones(1,nt); 
cc1=1;  [K1,M1]=fem.assema(gr,cc1,1,0); 
cc2=[1 0 0; 0 1 0; 0 0 1]; [K2,M1]=fem.assema(gr,cc2,1,0); 
max(max(abs(K1-K2))) 
%%
cc4=[ov1 0*ov1 0*ov1; 0*ov1 ov1 0*ov1; 0*ov1 0*ov1 1*ov1]; 
[K4,M1]=fem.assema(gr,cc4,1,0); max(max(abs(K1-K4)))
%%
cc5=[ov2 0*ov2 0*ov2; 0*ov2 ov2 0*ov2; 0*ov2 0*ov2 1*ov2]; 
[K5,M1]=fem.assema(gr,cc5,1,0); max(max(abs(K4-K5)))
%%
cc6=ov1; 
[K6,M1]=fem.assema(gr,cc6,1,0); max(max(abs(K6-K5)))