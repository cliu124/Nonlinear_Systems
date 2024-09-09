cmds1d:         script for 1D
bdmov1D:        creating a movie of the snaking branches 
cmds2dsq:       2D square domain 
cmds2dpsq:      perturbed square (still with D_4 symmetry)
cmds2dhex:      2D rectangular domain for hexagons 
cmds2dhexfro:   localized hex patterns on 2D long rectangular domain, 
cmds2dhexfroada: mesh adaptation of a solution from cmds2dhexfro,  
cmds2dhexb:     similar to cmds2dhex, but on domain twice as large; meant to 
                illustrate tips and tricks, for problems with ``too many'' patterns 
cmds2dsq_schemplots: 2D on square, to illustrate different bif-scenarios depending on nu 
bdmov2D:        creating a movie of the BD and sample solutions 
cmds2dtint:     patterns from initial guesses and time integration

cmds3dSC:       simple cube (SC) 3D patterns
cmds3dBCC:      body centered cube (BCC) 3D patterns 
cmdsBCClong:    localized BCC patterns from initial guesses
bdmov3D:        creating a movie of the localized 3D patterns 
cmdsBCClongref: mesh adaption for localized BCCs 


shinit, shinitpsq: init functions, with shinitpsq for the perturbed square
spjac:          Jacobian for fold continuation.

hfplot:         plot soln 1D with Hamiltonian
spl:            short-plot convenience function
spplots:        plot lattices and planforms (for tutorial) 
spplots3:       planforms in 3D (BCCs etc) 

etafu*:         different eta-functions for trullekrul mesh-adaption 