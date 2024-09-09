function opt=ilupparam(p,opt)
% ilupparam: checks entries in p.ilup and resets amgopt in lssAMG 
% opt=ilupparam(p,opt)
if isfield(p.ilup,'matching'); opt.matching=p.ilup.matching; end
if isfield(p.ilup,'ordering'); opt.ordering=p.ilup.ordering; end
if isfield(p.ilup,'droptol'); opt.droptol=p.ilup.droptol; end
if isfield(p.ilup,'droptolS'); opt.droptolS=p.ilup.droptolS; end
if isfield(p.ilup,'droptolc'); opt.droptolc=p.ilup.droptolc; end
if isfield(p.ilup,'condest'); opt.condest=p.ilup.condest; end
if isfield(p.ilup,'restol'); opt.restol=p.ilup.restol; end
if isfield(p.ilup,'maxit'); opt.maxit=p.ilup.maxit; end
if isfield(p.ilup,'elbow'); opt.elbow=p.ilup.elbow; end
if isfield(p.ilup,'lfil'); opt.lfil=p.ilup.lfil; end
if isfield(p.ilup,'lfilS'); opt.lfilS=p.ilup.lfilS; end
if isfield(p.ilup,'typetv'); opt.typetv=p.ilup.typetv; end
if isfield(p.ilup,'amg'); opt.amg=p.ilup.amg; end
if isfield(p.ilup,'npresmoothing'); opt.npresmoothing=p.ilup.npresmoothing; end
if isfield(p.ilup,'npostsmoothing'); opt.npostsmoothing=p.ilup.npostsmoothing; end
if isfield(p.ilup,'ncoarse'); opt.ncoarse=p.ilup.ncoarse; end
if isfield(p.ilup,'presmoother'); opt.presmoother=p.ilup.presmoother; end
if isfield(p.ilup,'postsmoother'); opt.postsmoother=p.ilup.postsmoother; end
if isfield(p.ilup,'FCpart'); opt.FCpart=p.ilup.FCpart; end
if isfield(p.ilup,'typecoarse'); opt.typecoarse=p.ilup.typecoarse; end
if isfield(p.ilup,'solver'); opt.solver=p.ilup.solver; end
if isfield(p.ilup,'damping'); opt.damping=p.ilup.damping; end
if isfield(p.ilup,'nrestart'); opt.nrestart=p.ilup.nrestart; end
if isfield(p.ilup,'mixedprecision'); opt.mixedprecision=p.ilup.mixedprecision; end
end