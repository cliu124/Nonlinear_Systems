% command templates for Gross-Pitaevski and vector GP; run cell-by-cell; 
close all; keep pphome; 
% first scalar GP 
p=[]; p=gpinit(p); p.nc.nsteps=30; %p=cont(p);
%% branch plotting 
plotbraf('p',3,3,'lab',[1,5]);
%% solution plotting 
splotsolf('p','pt1',4,1,5); pause; % VF plot
splotsolf('p','pt5',5,1,5); pause; % VF plot
splotsolf('p','pt5',6,10,3); pause; % amp plot
splotsolf('p','pt5',7,10,3); pause; % amp plot
splotsolf('p','pt5',8,13,2); % phase plot
%% vector GP 
close all; %clear all;
q=[];q=vgpinit(q); q=cont(q);
%% branch plotting 
figure(3);clf;
plotbra(q,3,3,'lab',[1,10,15]);
xlabel('\omega'); ylabel('p=max|u|/max|v| and p_1');
%% solution plotting
splotsolf('q','pt1',1,1,5); pause; % VF plot
splotsolf('q','pt10',1,1,5); % VF plot
%% look into phase invariance, i.e., check rotation, cond(G_u) and ker(G_u)
p=loadp('p','pt5');
splotsol(p,1,1,3);plotsol(p,2,2,3);
q=rotatesol(p,pi/2);
splotsol(q,3,1,3);plotsol(q,4,2,3);
zeroevcheck(p);
ocheck(p); % check orthogonality <Psi,G(U)>
%% some a posteriori error estimates
p=loadp('q','pt5');errcheck(p)