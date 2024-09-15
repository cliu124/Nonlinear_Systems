close all;keep pphome;
%% set parameters and run isc 
p=[]; p=ocinit(p,'c1','pt0','c1','pt0'); %p.tomopt.Order=6; 
p.u0(1:p.oc.s1.np)=1; % reset initial states to something 
p.u0(p.oc.s1.np+1:2*p.oc.s1.np)=1-p.oc.s1.u(p.oc.s1.nu+1); 
p.oc.T=10; p.oc.nti=31; alvin=0.25:0.25:1; p.oc.msw=1; p=isc(p,alvin);
%% cp and diagnostics plot
tit=['Eq to c1/pt0']; pcmp=1; v=[-30,30];
%tit=[fn.sd0 '/' fn.sp0 ' to ' fn.sd1 '/' fn.sp1];
psol3D(p.oc.s1,p.cp,4,pcmp,tit,0); view(v); colormap cool; 
switch pcmp; case 1; zlabel('v_1'); case 2; zlabel('v_2'); 
       case 3; zlabel('\lambda_1');  case 4; zlabel('\lambda_2'); end
zdia=lvdiagn(p,16,1); set(gca,'XTick',[0.1 1 10]); pause 
dia=lvdiagn(p,16,2); set(gca,'XTick',[0.1 1 10]); pause 
dia=lvdiagn(p,16,3); set(gca,'XTick',[0.1 1 10]); 