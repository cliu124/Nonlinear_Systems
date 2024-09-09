% script to test if ilupack is up 
A=[[1 2; 3 4]]; b=[1;2]; p=[]; 
try x=lssAMG(A,b,p); 
catch 
    fprintf('ilupack not installed/working: See the instructions in ' pphome '/README  \n'); 
end 