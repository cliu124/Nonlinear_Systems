function [p,cstop]=refufu(p,brout,ds)
% mod of STANUFU with test for refinement
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
eind=13; err=brout(6+eind); 
if err>p.nc.errbound; p.fuha.savefu(p); p=refineX(p,p.nc.sig); fprintf('\n');end 
cstop=0; p=retrigX(p); 
if(getlam(p)<p.nc.lammin)
    fprintf('  lam=%g < lammin=%g, stopping\n',getlam(p),p.nc.lammin); cstop=1;
end
if(getlam(p)>p.nc.lammax)
    fprintf('  lam=%g > lammax=%g, stopping\n',getlam(p),p.nc.lammax); cstop=1;
end
end
