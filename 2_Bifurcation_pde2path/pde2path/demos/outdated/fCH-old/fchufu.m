function cstop=fchufu(p,brout,ds)
% "user function" called after each cont.step 
brplot=brout(length(bradat(p))+p.plot.bpcmp); %y-axis value in bif-figure
lam=getlam(p);
if(p.sw.errchecksw>0) 
fprintf('%4i %6.4f  %6.4f  %5.3e  %5.3e  %4i  %s   %6.4f  ', ...
    p.file.count, lam, brplot, p.sol.res, p.sol.err, p.sol.iter, p.sol.meth, ds); 
else fprintf('%4i %6.4f  %6.4f  %5.3e  %4i  %s    %6.4f   ', ...
    p.file.count, lam, brplot, p.sol.res, p.sol.iter, p.sol.meth, ds); 
end 
for i=1:p.sw.npb; fprintf('%g  ',brout(i)); end; 
if(p.sw.spcalc==1) fprintf(' %2i ', p.sol.ineg); end; 
% put anything else here 
fprintf('\n');
cstop=0; 
if(lam<p.nc.lammin) 
    fprintf('lam=%g < lammin=%g, stopping\n',lam,p.nc.lammin); cstop=1; 
end 
if(lam>p.nc.lammax) 
    fprintf('lam=%g > lammax=%g, stopping\n',lam,p.nc.lammax); cstop=1; 
end 