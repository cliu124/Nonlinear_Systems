function stanheadfu(p)
% STANHEADFU: standard "header function" for runtime output called at start of cont/pmcont
%
%  stanheadfu(p)
%
% Output: "step lambda      y-axis    residual iter meth   ds"
%
% means: step # on branch, value of primary bif. param., y-axis value in bif diagram, 
%    residual, # newton iteration, continuation method (normal/arclength), stepsize
%
% See also cont, pmcont
auxdictl=0; if(isfield(p.plot,'auxdict')) auxdictl=length(p.plot.auxdict); end
if(auxdictl>=p.nc.ilam(1)) 
    primaux=p.plot.auxdict{p.nc.ilam(1)};
else primaux='lambda'; end
primaux(primaux=='\')=''; primaux=[primaux blanks(12-length(primaux))];
if(p.plot.bpcmp==0) yax='L2-norm';
elseif(p.plot.bpcmp<=auxdictl) yax=p.plot.auxdict{p.plot.bpcmp};
else yax='y-axis'; end
yax(yax=='\')=''; try; yax=[yax blanks(8-length(yax))]; catch; end; 
fprintf(char(['step   ' primaux yax 'residual  iter meth   ds       ']));
if(p.sw.errcheck>0) fprintf('err      '); end
if(p.sw.spcalc==1) fprintf('#-EV '); end 
npr=length(p.sw.bprint); 
for i=1:npr 
    if(p.sw.bprint(i)>0 && p.sw.bprint(i)<=auxdictl)
        addaux=char([p.plot.auxdict{p.sw.bprint(i)}]);
        addaux(addaux=='\')=''; fprintf(addaux);
        fprintf(blanks(11-length(addaux)));
    else fprintf('b(%i)       ',p.sw.bprint(i)); end
end
% put anything else here 
fprintf('\n');

