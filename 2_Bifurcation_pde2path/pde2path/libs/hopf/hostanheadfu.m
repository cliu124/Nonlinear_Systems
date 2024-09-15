function hostanheadfu(p)
% HOSTANHEADFU: standard "header function" for runtime output for Hopf 
auxdictl=0; if(isfield(p.plot,'auxdict')) auxdictl=length(p.plot.auxdict); end
if(auxdictl>=p.nc.ilam(1)) 
    primaux=p.plot.auxdict{p.nc.ilam(1)};
else primaux='lambda'; end
primaux(primaux=='\')=''; primaux=[primaux blanks(12-length(primaux))];
if(p.plot.bpcmp==0) yax='L2-norm';
elseif(p.plot.bpcmp<=auxdictl) yax=p.plot.auxdict{p.plot.bpcmp};
else yax='y-axis'; end
yax(yax=='\')=''; yax=[yax blanks(8-length(yax))];
fprintf(char(['step   ' primaux blanks(12-length(primaux)) yax ...
    blanks(11-length(yax)) ' res    iter  ds         T        nT   ']));
if(p.sw.errcheck>0) fprintf('err      '); end
if(p.sw.spcalc==1) fprintf('#-EV '); end 
npr=length(p.sw.bprint); addaux=''; 
for i=1:npr; 
    if(p.sw.bprint(i)>0 && p.sw.bprint(i)<=auxdictl)
        fprintf(char([p.plot.auxdict{p.sw.bprint(i)}]));
        addaux(addaux=='\')=''; fprintf(addaux);
        fprintf(blanks(11-length(primaux)));
    else fprintf('b(%i)       ',p.sw.bprint(i)); end
end
% put anything else here 
fprintf('\n');

