function [p,flag]=hoinistep(p,ds) 
% hoinistep: tomsol for first 2 points on branch; also comp.secant (para=3)
p0=p; 
% --------------- 1st point --------------------------------------
lam0=getlam(p); dlam=p.hopf.lam-lam0; % u0, u0dot set in hoswibra, here set new lam
p=setlam(p,p.hopf.lam); 
%p.hopf.tom.Stats='on'; p.hopf.tom.Stats_step='on';
[p,flag]=tomsol(p); 
if flag~=0; fprintf('hoinistep failed, try different ds\n'); p=p0; return; end 
fprintf('hoinistep got 1st point\n'); 
y0=p.hopf.y; t0=p.hopf.t; tl0=size(p.hopf.t); 
par=p.u(p.nu+1:end); f0=horhs(0,p.hopf.y(:,1),p,par); p.hopf.u0dot=f0(1:p.nu); 
p.sol.ptype=4; 
brout=[bradat(p); p.fuha.outfu(p,p.u)];          % userfu to append to bif-branches  
brplot=brout(length(bradat(p))+p.plot.bpcmp);    %y-axis value in bif-figure
p.branch=[p.branch brout];     
figure(p.plot.brfig); hold on; plot(getlam(p),real(brplot),'*');
[p,cstop]=p.fuha.ufu(p,brout,p.sol.ds); if(cstop==1); return; end  
p.hopf.tl=length(p.hopf.t); p.fuha.savefu(p); p.file.count=p.file.count+1; 
% ------------------------- 2nd point ----------------------------
lam0=getlam(p); lam1=lam0+dlam; p=setlam(p,lam1);  
p=tomsol(p); 
p.hopf.tl=length(p.hopf.t); 
if p.hopf.tl~=tl0; y0=interp1(t0,y0',p.hopf.t); y0=y0'; end % meshref in tom 
p.hopf.u0=p.hopf.y(1:p.nu,1); % save u and dot(u) for arclength 
f0=horhs(0,p.hopf.y(:,1),p,par); p.hopf.u0dot=f0(1:p.nu); 
p.hopf.ysec=[(p.hopf.y-y0); dlam*ones(1,length(p.hopf.t))]; 
p.hopf.ysec=p.hopf.ysec/honorm(p,p.hopf.ysec); % p.hopf.ysec(end-1:end,1), pause
p.hopf.y=[p.hopf.y; lam1*ones(1,length(p.hopf.t))]; % append lam-line 
brout=[bradat(p); p.fuha.outfu(p,p.u)];          % userfu to append to bif-branches  
brplot=brout(length(bradat(p))+p.plot.bpcmp);    %y-axis value in bif-figure
p.branch=[p.branch brout];   
figure(p.plot.brfig); hold on; plot(getlam(p),real(brplot),'*');
[p,cstop]=p.fuha.ufu(p,brout,p.sol.ds); if(cstop==1); return; end   
p.fuha.savefu(p); hoplot(p,p.plot.pfig,p.plot.pcmp,p.hopf.aux); 
p.file.count=p.file.count+1; fprintf('hoinistep done\n'); 