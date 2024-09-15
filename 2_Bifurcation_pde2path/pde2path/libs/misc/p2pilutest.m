%%
% p2pilutest: script to test if ilupack is up (and which one) 
A=sparse([1 2; 3 4]); b=[1;2]; p=stanparam; 
try; try; x=lssAMGb(A,b,p); fprintf('using ilupack from tu-BS\n'); 
    catch; x=lssILUb(A,b,p); fprintf('using ilupack from ilupack4m\n'); 
    end
catch 
    fprintf(['ilupack not installed/working: See the instructions in ' pphome '/README  \n']); 
end 