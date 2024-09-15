%% branch plot, b-cont, 
f=4; mclf(f); xlab='b'; c=8; ylab='E'; 
plotbra('b','pt20',f,c,'cl','k','lab',[1 20],'fms',0); 
plotbra('bb','pt25',f,c,'cl',p2pc('gr1'),'lab',[9 20],'fms',0);
xlabel(xlab); ylabel(ylab); grid on
%% soln plots 
global p2pglob; p2pglob.edc='none'; p2pglob.showbd=2; p2pglob.showN=0; 
pplot('b','pt1'); pause; pplot('b','pt20'); pause; pplot('bb','pt9'); pause; pplot('bb','pt20'); 
%% branch plot c0-cont, b=-1.66
f=3; mclf(f); xlab='c_0'; c=8; ylab='E'; %c=9; ylab='del_{mesh}'; 
plotbra('c00b','pt35',f,c,'cl','k','lab',[11 34]); 
plotbra('c00','pt15',f,c,'cl',p2pc('gr1'),'lab',14); 
xlabel(xlab); ylabel(ylab); grid on
%% soln plots
pplot('c00','pt14'); pause; pplot('c00','pt0'); pause; 
pplot('c00b','pt11'); pause; pplot('c00b','pt34');
%% branch plot c0-cont, b=-3.4 
f=4; mclf(f); xlab='c_0'; c=8; ylab='E'; %f=3; c=10; ylab='\int K dS'; %c=11; ylab='\delta_{mesh}'; 
plotbra('c01b','pt55',f,c,'cl','k','lab',55); %[15 30 55]); 
plotbra('c01','pt10',f,c,'cl',p2pc('gr1')); %,'lab',10); 
plotbra('c01b-1q','pt24',f,c,'cl','b','lab',[16 24]); 
plotbra('c01b-2','pt11',f,c,'cl','r','lab',9); 
plotbra('c01b-3','pt9',f,c,'cl','m','lab',16); 
xlabel(xlab); ylabel(ylab); grid on
%% soln plots 
p2pglob.showN=0; p2pglob.edc='none'; p2pglob.vi=[40,10]; 
pplot('c01b','pt55'); zticks('auto'); pause; pplot('c01b-2q','pt9'); pause; 
pplot('c01b-1q','pt16'); 
%%
p2pglob.edc='k'; p2pglob.vi=[40,20]; p2pglob.cm='parula'; 
pplot('c01b-1q','pt16'); pause; pplot('c01b-1q','pt24'); 
%% branch plot c0-cont, b=-4 
f=3; mclf(f); xlab='c_0'; c=8; ylab='E'; %c=9; ylab='del_{mesh}'; 
plotbra('c02b','pt32',f,c,'cl','k','fp',7); 
plotbra('c02b-1q','pt20',f,8,'cl','b','lab',19); 
plotbra('c02b-2q','pt18',f,c,'cl','r','lab',18); 
plotbra('c02b-2q-1','pt10',f,c,'cl','m','lab',10); 
xlabel(xlab); ylabel(ylab); grid on
%% soln plots 
global p2pglob; p2pglob.edc='none'; p2pglob.cut=0; p2pglob.showbd=2; 
p2pglob.showN=0; p2pglob.cb=0; p2pglob.axlab=0; p2pglob.vi=[20,10]; 
pplot('c02b-1q','pt19'); pause; pplot('c02b-2q','pt18');pause; pplot('c02b-2q-1','pt10');
%% results for b=-4; 
f=3; mclf(f); xlab='c_0'; c=8; ylab='E'; %c=9; ylab='del_{mesh}'; 
plotbra('c02b','pt32',f,c,'cl','k','fp',7); 
plotbra('c02b-1q','pt28',f,8,'cl','b'); %,'lab',26); 
plotbra('c02b-2q','pt30',f,c,'cl','r','lab',30); 
%plotbra('rt','pt20',f,c,'cl',p2pc('g1'),'lab',28); 
plotbra('c02b-2q-1','pt20',f,c,'cl',p2pc('g1'),'lab',[10,20]); 
xlabel(xlab); ylabel(ylab); grid on
%%
global p2pglob; p2pglob.edc='none'; p2pglob.cut=0; p2pglob.showbd=2; 
p2pglob.showN=0; p2pglob.cb=0; p2pglob.axlab=0; p2pglob.vi=[5,40]; 
pplot('c02b-1q','pt26'); pause; pplot('c02b-2q','pt30');pause; pplot('c02b-2q-1','pt10');
pause; pplot('c02b-2q-1','pt20');
%%
p2pglob.edc='k';p2pglob.vi=[20,5]; %pplot('1/c0b-1q','pt10'); pause;
pplot('1/c0b-1q','pt15'); pause; 
pplot('1/c0b-2q','pt16'); pause; pplot('1/c0b-2q-1','pt12'); 
%%
p=loadp('1/c0b-2q','pt16'); pplot(p); out=hcapbra(p,p.u); pause
p=loadp('1/c0b','pt60'); pplot(p); out=hcapbra(p,p.u); pause 
p=loadp('1/c0','pt15'); pplot(p); out=hcapbra(p,p.u); 
%%
f=3; mclf(f); xlab='c_0'; c=10; ylab='E'; %c=9; ylab='del_{mesh}'; 
plotbra('1/c1b','pt60',f,c,'cl','k','lab',[20 50]); 
plotbra('1/c1','pt40',f,c,'cl','b','lab',[10 15]); 
plotbra('1/c1b-1','pt20',f,c,'cl','r','lab',[10 15]); 
%%
p=loadp('1/c1b','pt20'); pplot(p); out=hcapbra(p,p.u); pause
p=loadp('1/c0b','pt50'); pplot(p); out=hcapbra(p,p.u); pause 
p=loadp('1/c0','pt15'); pplot(p); out=hcapbra(p,p.u); 

