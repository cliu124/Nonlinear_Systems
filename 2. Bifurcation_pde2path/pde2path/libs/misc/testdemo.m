function testdemo(ddir,cmds)
% testdemo(ddir,cmds): run ddir/cmds.m and write log.txt to calling dir 
global pphome; 
cd(pphome); fID=fopen('log.txt','a'); fprintf(fID,['testing ' ddir '/' cmds ' ... ']); 
warning('off','all'); 
cd([pphome '/libs/misc']); copyfile nopause.m pause.m; cd([pphome ddir]); 
try run(cmds); cd(pphome); fID=fopen('log.txt','a'); fprintf(fID,'worked\n'); 
catch; cd(pphome); fID=fopen('log.txt','a'); fprintf(fID,'failed\n'); 
end
cd([pphome '/libs/misc']); delete pause.m; cd(pphome); warning('on','all'); 
end