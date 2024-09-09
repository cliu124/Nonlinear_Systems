% script pubidx
% : publish library index html files (internal use)
pubopt=[]; pubopt.format = 'html'; pubopt.showCode = true;
pubopt.createThumbnail=false; pubopt.useNewFigure=false;
pubopt.evalCode=false;
pubopt.outputDir='/home/hu/path/pde2path/html/index'; 
pubopt.stylesheet=''; %'/home/hu/path/internal/mystyle.xsl';
sdir=pwd; 
cd('/hh/path/internal/index'); 
aa=what(pwd); totalfiles=numel(aa.m), fileName=aa.m; 
for i=1:1:totalfiles; publish(fileName{i},pubopt); end
cd('/hh/path/pde2path/libs/p2p/'); publish('stanparam.m',pubopt);
cd('/hh/path/pde2path/libs/hopf/'); publish('hostanparam.m',pubopt);
cd(sdir)