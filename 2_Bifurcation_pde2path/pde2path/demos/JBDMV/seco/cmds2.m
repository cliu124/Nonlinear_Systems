%% spatially homogeneous primary Hopf; use tl=40 points in t; 
ds=0.1; aux.tl=40; p=hoswibra('0','hpt1',ds,4,'H1',aux); p.nc.dsmax=2;
p.sw.bifcheck=0; p.file.smod=1;% switch off bif.detection, save every point 
p.hopf.flcheck=1; p.hopf.fltol=1e-6; % switch on multiplier comp  
p.sw.verb=0; p=setbel(p,2,1e-4,5,@lss); p=cont(p,40); % use bel, and go! 
%% splice together Hopf and Turing: choose xcut and j0 by trial and error
p1=loadp('T1','pt25'); p=loadp('H1','pt25'); xcut=-160; j0guess=3.35; 
x=getpte(p);idx=find(x>xcut);%find indizes where to replace Hopf by Turing
p.u(p.nu+1)=j0guess; p.hopf.lam=j0guess; % set j0 to the guess j0guess 
for j=1:p.hopf.tl % put Turing-soln into Hopf-data 
 p.hopf.y(idx,j)=p1.u(idx); p.hopf.y(idx+p.np,j)=p1.u(idx+p.np); 
end 
for i=1:1;  hoplot(p,1,i); end; pause % plot for checking 
p.sw.verb=2; p.sol.ds=-0.01; p.nc.dsmin=0.01; 
p.nc.imax=20; p.nc.tol=1; p.hopf.flcheck=0; p.nc.dsmax=5; 
p=resetc(p); p=setfn(p,'lh1a'); p.hopf.nfloq=10; 
p=cont(p,1); % 1 step to see if initial guess was OK 
%% a few more steps to get on the branch 
p=resetc(p); p.nc.imax=10; p.nc.tol=1e-2; p=cont(p,3); p.nc.tol=1e-4; p=cont(p,3); 
%% continue ! 
p=resetc(p); 
p.nc.imax=10; p.nc.tol=1e-6; p.hopf.flcheck=0; p=cont(p,1); p.file.smod=10; p=cont(p,50); 
%% other direction 
p=loadp('lh1a','pt1','lh1b'); p.sol.ds=-p.sol.ds; p=resetc(p); p=cont(p,200); 
%% BD 
fnr=3; figure(fnr); clf; cmp=8; % desired branch-compo,  
ptli0=[16 22]; ptli1=[1 40]; ptli2=[20 90 110 170 230 290]; % label lists 
plotbra('lh1a','pt40',fnr,cmp,'cl','r','lab',ptli1,'lp',40); 
plotbra('lh1b','pt300',fnr,cmp,'cl','m','lab',ptli2,'fp',4); 
plotbra('H1','pt50',fnr,cmp,'cl','b','lab',ptli0,'fp',16); 
xlabel('j0'); ylabel('||U||'); 
%% soln plots, with Fl-multipl., uncomment for desired directory 
v=[15 60]; aux.nfloq=30;  % view, and number of Floquet-multipliers 
dir='H1';ptli=ptli0; 
%dir='lh1a';ptli=ptli1; 
%dir='lh1b';ptli=ptli2; 
for i=ptli; 
    pt=['pt' mat2str(i)]; hoplotf(dir,pt,1,1); xlabel('t'); 
    figure(1); title([dir '/' pt]); ylabel('t'); zticks([4 8]); 
    view(v); shading interp;  muv1=floqap(dir,pt,aux);  
    pause
end