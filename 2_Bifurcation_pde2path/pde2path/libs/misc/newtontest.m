function p=newtontest(dir,pt,varargin)
% newtontest: test different Newton loop settings: standard, chord, fsolve, nleq1
del=0.1; xtol=1e-3; fs=1; % do fsolve by default 
nleq1sw=0; % don't do nleq1 by default
if nargin>2; del=varargin{1}; end 
if nargin>3; xtol=varargin{2}; end 
if nargin>4; fs=varargin{3}; end 
if nargin>5; nleq1sw=varargin{4}; end 
p=loadp(dir,pt); p.u(1:p.nu)=p.u(1:p.nu)+del; plotsol(p); r=resi(p,p.u);
ftol=p.nc.tol; 
fprintf('ini-res=%g\n',norm(r,'inf')); p0=p; 
p.sw.newt=0; p.nc.tol=1e-8; % standard Newton 
tic; [u1,r1,i1,Gu,Glam,p]=nloop(p,p.u); t1=toc; 
fprintf('standard newton: res=%g, iter=%g, time=%g\n',r1,i1,t1); 
p.u=u1; plotsol(p,2,1,p.plot.pstyle); 
p=p0; p.sw.newt=0; p.nc.tol=1e-8; % chord 
tic; [u1b,r1b,i1b,Gu,Glam,p]=nloop(p,p.u); t1b=toc; 
fprintf('chord: res=%g, iter=%g, time=%g\n',r1b,i1b,t1b); 
p.u=u1; plotsol(p,4,1,p.plot.pstyle); 
if fs>0; p=p0; p.fsol.fsol=1; p.nc.tol=ftol; p.fsol.imax=10; 
tic; [u3,r3,i3,Gu,Glam,p]=nloop(p,p.u); t3=toc; 
fprintf('fsolve: res=%g, iter=%g, time=%g\n',r3,i3,t3); 
p.u=u3; plotsol(p,5,1,p.plot.pstyle); 
end 
if nleq1sw>0; p=p0; p.sw.newt=2; p.fsol.fsol=0; p.nc.tol=xtol; % NLEQ1 
tic; [u2,r2,i2,Gu,Glam,p]=nloop(p,p.u); t2=toc; 
fprintf('NLEQ1: res=%g, iter=%g, time=%g\n',r2,i2,t2); 
p.u=u2; plotsol(p,6,1,p.plot.pstyle); 
end  
fprintf('standard newton (Fig.2): res=%g, iter=%g, time=%g\n',r1,i1,t1); 
fprintf('chord (Fig.4): res=%g, iter=%g, time=%g\n',r1b,i1b,t1b); 
if fs>0; fprintf('fsolve (Fig.5): res=%g, iter=%g, time=%g\n',r3,i3,t3); end 
if nleq1sw>0; fprintf('NLEQ1 (Fig.6): res=%g, iter=%g, time=%g\n',r2,i2,t2); end 