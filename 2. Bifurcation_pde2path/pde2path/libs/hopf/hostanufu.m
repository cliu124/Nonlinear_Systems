function [p,cstop]=hostanufu(p,brout,ds)
% STANUFU: standard "user function" called after each cont.step, HOPF 
p=ulamcheck(p); % test if a desired lambda value has been passed
brplot=brout(length(bradat(p))+p.plot.bpcmp); %y-axis value in bif-figure
T=0; 
switch p.sw.para;  % period T treated differently depending on para
    case 3; % TOM 
        T=p.hopf.T; 
try Td=p.hopf.ysec(2*p.nu+1,1); lamd=p.hopf.ysec(p.nu+2,1); catch Td=0; lamd=0; end;
    case 4; % TOM for assembling 
        T=p.hopf.T; 
try; Td=p.hopf.tau(end-1); lamd=p.hopf.tau(p.nc.neq*p.np*p.hopf.tl+2); 
catch Td=0; lamd=0; end;
    case 6; % NKrylov 
        T=p.u(p.nu+p.hopf.iT); 
try Td=p.hopf.tau(end-1); lamd=p.hopf.tau(p.nu+2); catch Td=0; lamd=0; end; 
    otherwise; Td=0; lamd=0; 
end 
fprintf('%4i %s %s %s %2i %s %s %4i %s %s', ...
    p.file.count, printcon(getlam(p)), printcon(brplot), printcon(p.sol.res), p.sol.iter, ...
    printcon(ds), printcon(T), length(p.hopf.t), printcon(Td), printcon(lamd));
npr=length(p.sw.bprint);
for i=1:npr; fprintf('%s ',printcon(brout(length(bradat(p))+p.sw.bprint(i)))); end;
% put anything else here
fprintf('\n'); cstop=0;
if(getlam(p)<p.nc.lammin)
    fprintf('  lam=%g < lammin=%g, stopping\n',getlam(p),p.nc.lammin); cstop=1;
end
if(getlam(p)>p.nc.lammax)
    fprintf('  lam=%g > lammax=%g, stopping\n',getlam(p),p.nc.lammax); cstop=1;
end
end
