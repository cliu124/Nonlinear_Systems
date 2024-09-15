function [p,cstop]=coarsufu(p,brout,ds)
% mod of STANUFU with coarsening 
p=ulamcheck(p); % test if a desired lambda value has been passed
brplot=brout(length(bradat(p))+p.plot.bpcmp); %y-axis value in bif-figure
fprintf('%4i %s %s %5.2e %4i  %s %s ', ...
    p.file.count,printcon(getlam(p)),printcon(brplot),p.sol.res,p.sol.iter, ...
    p.sol.meth,printcon(ds));
if(p.sw.errcheck>0) fprintf('%5.2e ',p.sol.err);end
if(p.sw.spcalc==1) fprintf(' %2i ', p.sol.ineg); end;
npr=length(p.sw.bprint); 
for i=1:npr; fprintf('%s ',printcon(brout(length(bradat(p))+p.sw.bprint(i)))); end;
fprintf('\n');
% put anything else here
delmind=11; delm=brout(6+delmind); keepbd=1; sig=p.nc.sig; nit=6; 
if delm>p.nc.errbound; p.fuha.savefu(p); p=degcoarsenX(p,sig,nit,keepbd); 
   p.sol.ds=0.01;  fprintf('\n');end 

cstop=0;
if(getlam(p)<p.nc.lammin)
    fprintf('  lam=%g < lammin=%g, stopping\n',getlam(p),p.nc.lammin); cstop=1;
end
if(getlam(p)>p.nc.lammax)
    fprintf('  lam=%g > lammax=%g, stopping\n',getlam(p),p.nc.lammax); cstop=1;
end
end
