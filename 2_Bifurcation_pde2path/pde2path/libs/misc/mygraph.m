function G=mygraph(np,sw)
% mygraph: convenience function to generate different graphs (or load from disk) 
% 
% sw=0: load from disk, np is filename   (else np is #of nodes in graph) 
%    1: 'hand-made' ((just for testing) 
%  2,3: Watts-Strogatz, short average path lengths and high clustering
%    4:  Barabasi–Albert (BA) model. Scale free, little soln clustering 
switch sw; 
    case 0; G=load(np);  G=G.G; % load Graph from disk 
    case 1; s=[1 1 1 1 1 1 1 2 3 4 5 6 7 8]; 
            t=[2 3 4 5 6 7 8 3 4 5 6 7 8 2]; G=graph(s,t); 
    case 2; K=4; beta=0.2; G=WattsStrogatz(np,K,beta); % WS 
    case 3; K=2; beta=0.8; G=WattsStrogatz(np,K,beta); 
    case 4; mm=3; A=BAgraphA(np, mm, mm); G=graph(A);  % Barabasi–Albert        
end
h=plot(G); h.NodeLabel={}; 