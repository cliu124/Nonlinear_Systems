function plotall1(p) % convenience func for plotting abs/real/imag/VF, calling userplot
p.plot.pcmp=10; plotsol(p,1,10,-1); pause; % abs, real, imag 
p.plot.pcmp=13; plotsol(p,1,10,-1); figure(1); view(0,90); pause % angle 
p.plot.pcmp=14; plotsol(p,1,10,-1); % quiver 