% biocaps: geomtut, §4.2; Helfrich SC model with stress free BCs
cmds1.m:       c0=0.5; al=1; b=-0.5 (initial), then cont in c0 at different b
cmds1plot.m:   plotting for cmds1 results
cmds2.m:       cont of non-axi branches from cmds1 in b and in lam1 

hcapbra.m:     [pars; V; A; E1; E2; bdE1; meshqual]
%               1-5   6  7  8   9   10   11
bdint.m:       boundary integral, trapezoidal rule  