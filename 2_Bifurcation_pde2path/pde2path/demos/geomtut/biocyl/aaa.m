% Helfrich cylinder, with clampded BCs, geomtut appendix B
cmds1.m:   cont in c0
cmds2.m:   cont in al
cmds3.m:   cont in c0 at al=1.4 (b/*) and at al=0.6 (s/*)
cmds4.m:   cont in lam1 (c0=0.7)
cmds4b.m:  cont in lam1 (c0=0)
hcylbra.m: branch output
hcylbra2.m: appending lam1*A to the branch output (afterthought, see 
   cmdsaux.m for how to step through branches and store and plot the
   additional output) 
cylinit.m: init


bdmov*.m:  scripts to do movies of BDs 