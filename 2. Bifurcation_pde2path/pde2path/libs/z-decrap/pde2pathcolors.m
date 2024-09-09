function v=pde2pathcolors(str)
% pde2pathcolors: convenience function to choose colors. 
%
% v=pde2pathcolors(str)
% Available colors include 
% g* (green), b* (blue), y* (yellow), r* (red), v* (violet)
% br* (brown), o* (orange), gr* (grey),  where *=1,2,3 
switch str
    case 'g1'; v=[0 100 0]/255;
    case 'g2'; v=[110,139,61]/255;
    case 'g3'; v=[152,251,152]/255;    
    
    case 'b1'; v=[0,0,139]/255;
    case 'b2'; v=[72,61,139]/255;
    case 'b3'; v=[16,78,139]/255;
        
    case 'y1'; v=[255,185,15]/255;
    case 'y2'; v=[255,215,0]/255;
    case 'y3'; v=[238,221,130]/255;
    
    case 'r1'; v=[139,0,0]/255;
    case 'r2'; v=[205,79,57]/255;
    case 'r3'; v=[255,48,48]/255;
        
    case 'v1'; v=[139,0,139]/255;
    case 'v2'; v=[205,0,205]/255;
    case 'v3'; v=[125,38,205]/255;
    
    
    case 'br1'; v=[139,69,19]/255;
    case 'br2'; v=[210,105,30]/255;
    case 'br3'; v=[139,90,43]/255;
        
    case 'o1'; v=[255,140,0]/255;
    case 'o2'; v=[255,165,0]/255;
    case 'o3'; v=[238,121,66]/255;
    
    
    case 'gr1'; v=[105,105,105]/255;
    case 'gr2'; v=[119,136,153]/255; 
    case 'gr3'; v=[168,168,168]/255;
        
    case 'yellow'; v=[1 1 0];
    case 'magenta'; v=[1 0 1];
    case 'cyan'; v=[0 1 1];
    case 'red'; v=[1 0 0];
    case 'green'; v=[0 1 0];
    case 'blue'; v=[0 0 1];
    case 'white'; v=[1 1 1];
    case 'black'; v=[0 0 0];
        
        
end
          
 % if matlab standard colors like 'b' or 'r' are used       
 if length(str)==1;v=bitget(find('krgybmcw'==str)-1,1:3);end
    
