% demoscript for Gross-Pitaevski; 
close all; clear all; 
fprintf('gpinit: first define dummy starting guess on regular mesh,\n');
fprintf('use meshref to generate fine mesh in middle of domain, \n');
fprintf('and finally create starting guess on the finer mesh and use \n');
mypause('meshref again to create starting point on suitable mesh');
p=[];p=gpinit(p);
mypause('now run continuation'); 
p=cont(p);
% ---------------------------------------------------------
fprintf('Plot branch and some solutions in different plotstyles\n'); 
mypause('see plotsol.m for the customized plotting');
plotbra(p,3,3,'lab',[10,15]);
plotsolf('p','p1',4,1,5); % VF plot
plotsolf('p','p10',5,1,5); % VF plot
plotsolf('p','p10',6,10,3); % amp plot
plotsolf('p','p15',7,10,3); % amp plot
plotsolf('p','p20',8,13,2); % phase plot
% --------------------------------------------------------- 
mypause('vgpinit: proceed similar to gpinit to get starting point');
q=[];q=vgpinit(q);
mypause('Continuation with mesh-adaption every 3 steps:'); 
q.amod=3;q=cont(q);
% --------------------------------------------------------- 
mypause('plot p and q branch together, and some q solutions'); 
clf(3);
plotbra(p,3,3,'lab',[1,10,15]);
plotbra(q,3,3,'lab',[1,10,15],'cl','r');
xlabel('\omega'); ylabel('p=max|u|/max|v| and p_1');
plotsolf('q','p1',10,1,5); % VF plot
plotsolf('q','p10',11,1,5); % VF plot