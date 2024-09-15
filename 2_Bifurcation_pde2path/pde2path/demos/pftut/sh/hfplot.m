function hfplot(dir,pt) % convenience function for plotting soln, H, and Fourier
p=loadp(dir,pt); clf(1); geth(p,p.u); fouplot(p,10,1,25,0,[dir '/' pt]); 
