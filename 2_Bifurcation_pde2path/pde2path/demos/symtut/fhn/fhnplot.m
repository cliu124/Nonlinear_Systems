%%
figure(3); clf; plotbra('prep/decoup','pt20',3,0,'lab',[10,20]);
%%
plotsol('prep/decoup','pt10'); pause; plotsol('prep/decoup','pt20'); pause; 
plotsol('prep/cubic','pt5'); 
%% cell 8 - plot the supercritical pitchfork and sols
figure(3); clf; plotbra('stand','pt8',3,2,'lab',7);
plotbra('travel','pt11',3,2,'cl','b','lab',11);
plotbra('travel-','pt11',3,2,'cl','m','lab',11);
axis([0.7 1.5 -0.02 0.02]); 
%% 
clf(3); plotbra('asym',3,2,'lab',2); 
%%
clf(3); plotbra('cusp',3,3,'labi',9); 
%%
plotsol('stand','pt7'); axis([-10 10 -1.1 1.1]); pause; 
plotsol('travel','pt11'); axis([-10 10 -1.1 1.1]);