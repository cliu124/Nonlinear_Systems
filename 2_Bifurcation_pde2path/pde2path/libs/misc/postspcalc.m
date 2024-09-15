function il=postspcalc(varargin)
% POSTSPCALC: to be used when p.sw.spcalc was 0 and one a posteriori wants 
% to compute stability. Works only if p.file.smod=1, i.e. every solution is
% saved.
%
% il=postspcalc('dir'): dir is the folder of the solutions for which one 
%                       wants to calc. the stability afterwards. 
%                       The stablity information can be found in il. 
%                       Furthermore, it is stored in p.branch of the last 
%                       point, which is saved in the folder dir.
%
% il=postspcalc('dir',aux): aux must be a structure. Possible fields are
%
% aux.neig:     replaces p.nc.neig by aux.neig.
% aux.num:      calculates the number of positive eigenvalues for labels, 
%               which are stored in aux.num. aux.num must be a row vector.
%               For instance num=[10:20].
if nargin==1
    aux=[];
end
if nargin==2
    fol=varargin{1}; %source folder 
    aux=varargin{2};
end
    fol=varargin{1}; %source folder 
if isfield(aux,'num')
    num=aux.num;
else
    num=sort(getlabs(fol))'; % labels of file in the source folder
end
il=-1*ones(1,length(num)); %ineg line in p.branch
parfor j=1:length(num)
    i=num(j); pp=loadp(fol,['pt',num2str(i)]);
    if isfield(aux,'neig');pp.nc.neig=aux.neig;end
    r=resi(pp,pp.u); Gu=getGu(pp,pp.u,r); [ineg,~]=vspcalc(Gu,pp);
    il(j)=ineg; fprintf('j=%i, ineg=%i\n',i,ineg); 
end
p=loadp(fol,['pt',num2str(num(end))]);
sil=numel(il);
p.branch(3,end-sil+1:end)=il;
p.file.count=p.file.count-1;
p.fuha.savefu(p);