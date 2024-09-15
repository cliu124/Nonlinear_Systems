%% branch plot
f=3; c=1; mclf(f); %plotbra('sNb',f,c,'cl',p2pc('gr2'),'lab',4);
plotbra('N','pt60',f,c,'cl','k','lab',52);
plotbra('Nb','pt2',f,c,'cl',p2pc('gr1'),'lab',2);
plotbra('Nr1','pt12',f,c,'cl',p2pc('gr1'),'lab',[2 12]);
plotbra('Nr2','pt10',f,c,'cl',p2pc('gr2'),'lab',[]);
plotbra('Nr3','pt15',f,c,'cl',p2pc('gr1'),'lab',[15]);;
plotbra('N1','pt29',f,c,'cl','b','lab',20); 
plotbra('N2',f,c,'cl',p2pc('r'),'lsw',0,'lab',25); 
plotbra('N3','pt15',f,c,'cl',p2pc('r2'),'lab',10); 
plotbra('N3-1',f,c,'cl','m','lsw',0,'lab',[20]); 
plotbra('N4',f,c,'cl',p2pc('o1'),'lab',15); 
plotbra('N5',f,c,'cl',p2pc('g1'),'lab',8);
plotbra('N6','pt8',f,c,'cl',p2pc('r3'),'lab',6);
%% branch plot, zoom 
f=4; c=1; mclf(f); plotbra('N','pt35',f,c,'cl','k','lab',[],'fp',10);
plotbra('N2',f,c,'cl','r','lab',222,'fp',38); 
plotbra('N3',f,c,'cl',p2pc('r1'),'lab',13); 
plotbra('N3-1',f,c,'cl','m','lab',20); 
plotbra('N4',f,c,'cl',p2pc('o1'),'lab',16); 
plotbra('N5',f,c,'cl',p2pc('g1'),'lab',8);
axis([41 75 0.55 0.74]); xlabel(''); ylabel(''); 
%% soln plot
p2pglob.tsw=2; p2pglob.vi=[20 20]; p2pglob.cut=6; p2pglob.edc='k'; p2pglob.showbd=2; 
plotsol('N','pt52'); pause; plotsol('Nr1','pt2'); pause; 
plotsol('Nr1','pt12'); pause; plotsol('Nr2','pt10'); pause; 
plotsol('Nr3','pt15'); 
%%
p2pglob.edc='none'; p2pglob.cut=1; p2pglob.cb=0; 
plotsol('Nb','pt2'); pause; plotsol('N','bpt1'); pause; 
plotsol('Nr3','pt15'); pause; 
plotsol('N1','pt20'); pause; plotsol('N2','pt25'); pause; 
plotsol('N3','pt10'); pause; plotsol('N3-1','pt20'); pause; 
plotsol('N4','pt15'); pause; plotsol('N5','pt8'); pause; plotsol('N6','pt6');
%%
p2pglob.edc='k';  p2pglob.cb=0; p2pglob.cut=0; plotsol('N1','pt6'); 