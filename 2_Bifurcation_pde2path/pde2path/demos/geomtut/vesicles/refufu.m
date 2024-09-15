function [p,cstop]=refufu(p,brout,ds)
global p2pglob;
% mod of STANUFU with refinement (by area AND by triangle shape) and degcoarsenX
p=ulamcheck(p); % test if a desired lambda value has been passed
brplot=brout(length(bradat(p))+p.plot.bpcmp); %y-axis value in bif-figure
fprintf('%4i %s %s %5.2e %4i  %s %s ', ...
    p.file.count,printcon(getlam(p)),printcon(brplot),p.sol.res,p.sol.iter, ...
    p.sol.meth,printcon(ds));
if(p.sw.errcheck>0) fprintf('%5.2e ',p.sol.err); end
if(p.sw.spcalc==1) fprintf(' %2i ', p.sol.ineg); end
npr=length(p.sw.bprint); 
for i=1:npr; fprintf('%s ',printcon(brout(length(bradat(p))+p.sw.bprint(i)))); end
fprintf('\n');
% put anything else here
try npmax=p2pglob.npmax; catch; npmax=3500; end % check for (soft) upper bound on np 
p.ref=0; % switch used in sscontrol 
mq=meshqdat(p); err=mq(1); count=0; res=inf;
while err>p.nc.delbound && count<10 && res>p.nc.tol; % degcoarsen-refine-retrig
    p.fuha.e2rs=@e2rsshape1; 
    nt=p.nt; p=degcoarsenX(p,p.nc.sigc,100); 
    [u1,res,~]=nloop(p,p.u); up=p.up;
    [p,u1]=updX(p,u1); p.u=u1; p.up=up;
    if p.np>npmax;  sig=nt/p.nt-1; sigr=sig/3; p=refineX(p,sigr); 
    else, p=refineX(p,p.nc.sigr); end; 
    p=retrigX(p); 
    [u1,res,~]=nloop(p,p.u); up=p.up;
    [p,u1]=updX(p,u1); p.u=u1; p.up=up;
    p.sol.xi=1/p.nu; mq=meshqdat(p); err=mq(1);
 p.ref=1; count=count+1;
end 
A=doublearea(p.X,p.tri); Ai=sort(A/2,'descend'); p.nc.rlong=1; 
if Ai(1)>p.nc.Ab; p.fuha.e2rs=@e2rsA; s=find(Ai>p.nc.Ab); % refine large area  triangles 
    try; sf=p.nc.sigfac; catch sf=1; end; sig=sf*length(s)/p.nt; 
  p=refineX(p,sig); p=retrigX(p); 
  [u1,res,~]=nloop(p,p.u); up=p.up;
    [p,u1]=updX(p,u1); p.u=u1; p.up=up;
  p.sol.xi=1/p.nu; 
  p.ref=1; 
end
cstop=0; % further proceed as usual 
if(getlam(p)<p.nc.lammin)
    fprintf('  lam=%g < lammin=%g, stopping\n',getlam(p),p.nc.lammin); cstop=1;
end
if(getlam(p)>p.nc.lammax)
    fprintf('  lam=%g > lammax=%g, stopping\n',getlam(p),p.nc.lammax); cstop=1;
end