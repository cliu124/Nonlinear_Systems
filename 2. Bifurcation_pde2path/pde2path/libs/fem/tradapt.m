function [tri,xy,Nmetric,bndmesh,triQ]=tradapt(tri,xy,Nmetric,bndmesh,geomfunc,options)
% tradapt: minor mod of trullekrul adapt_mesh with switch sw=options.sw
% binary coding, 
%       sw=scrm for swap/coarsen/refine/move, 
% i.e., m=0/1 -> moving off/on, r=0/1 -> refine off/on etc 
% hence, decimal: 
%  sw       1      2      3      4       5    6     7     8     9  ... 15
% action: move   refine  m,r  coarsen   c,m  c,r  c,r,m  swap  s,m    all  
% 
% hence sw=15: orginal adapt_mesh behaviour
%
% typically, adaptm is called by oomeshada and oomeshadactr, and 
% options=p.trop   (trullekrul-options) 
%
if (size(xy,2) == 2 && size(Nmetric,2) == 3) || (size(xy,2) == 3 && size(Nmetric,2) == 6)
triQ=elem_qual(tri, xy, Nmetric, options);
elseif size(xy,2)
triQ=elem_qual(tri, xy, Nmetric(:,1:3), options);
else
triQ=elem_qual(tri, xy, Nmetric(:,1:6), options);	
end;

minTriQ=min(triQ);
if minTriQ < options.minqual
	warning(sprintf('Low quality elements in input(%0.0e)',minTriQ));
	options.minqual=minTriQ;
	if minTriQ < 0
		error('inverted element in input mesh');
	end;
end;

if isfield(bndmesh,'edg') || isfield(bndmesh,'fac')
	if isfield(bndmesh,'crnds')
	bndmesh.crnds=[bndmesh.crnds; geom_crnds(bndmesh,1)];
	else
	bndmesh.crnds=geom_crnds(bndmesh,1);
	end; 
	bndmesh.crnds=unique(bndmesh.crnds);
end;

if size(xy,2) == 2 && not(isfield(bndmesh,'edg'))
	bndmesh=bndmesh_polygon(tri,xy,bndmesh,options);
elseif size(xy,2) == 3 && not(isfield(bndmesh,'fac'))
	bndmesh=bndmesh_polyhedron(tri,xy,bndmesh,options);
end;
if options.mntn_bks
	bks=bks_init(tri);
else
	bks=[];
end;
if options.debug
     sanity_check(tri,xy,triQ,Nmetric,options);
end;
%options.prag_adapt, pause 
if options.prag_adapt == 3 %only smooth
[xy,triQ,bks,ndone4]=adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options,options.innerit);
elseif options.prag_adapt
[tri,xy,bndmesh,Nmetric,triQ,bks] =pragadapt(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
else   
[tri,xy,bndmesh,Nmetric,triQ,bks] =newadapt(tri,xy,Nmetric,bndmesh,triQ,bks,bks,geomfunc,options);
end;
if options.debug
     sanity_check(tri,xy,triQ,Nmetric,options);
end;

function [tri,xy,bndmesh,Nmetric,triQ,bks] =pragadapt(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options)	
  save_plot(tri,xy,0,options,'initial mesh');  ndone1=0; ndone2=0; ndone3=0; ndone4=0; 
  sw=options.sw; % HU, switch encoding behaviour 
  for i=1:options.innerit;  if options.verbose == 1;  disp_qual(triQ,tri,xy,i); end;
    if ismember(sw,[4:7, 12:15]) %COARSENING
        if options.RMedg            
        [xy,tri,bndmesh,Nmetric,triQ,bks,ndone1]=adapt_rm_edg(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
        else
      	[tri,triQ,bks,bndmesh,ndone1,xy,Nmetric]=adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options);
        end;
        save_plot(tri,xy,(i-1)*4+2,options,sprintf('i=%0.0f, after coarsening',i));
    end
    if ismember(sw,8:15) %SWAPPING
        [tri,triQ,bks,bndmesh,ndone2]=adapt_flipedg(tri,xy,Nmetric,bndmesh,triQ,bks,options); 
        save_plot(tri,xy,(i-1)*4+1,options,sprintf('i=%0.0f, after swapping',i));
    end
	if ismember(sw,[2,3,6,7,10,11,14,15]); %REFINEMENT
        [xy,tri,bndmesh,Nmetric,triQ,bks,ndone3]=adapt_add_nd(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
        save_plot(tri,xy,(i-1)*4+3,options,sprintf('i=%0.0f, after refinement',i));
    end
	if ismember(sw,[1,3,5,7,9,11,13,15]); %SMOOTHING
       if options.prag_adapt == 2 || (options.prag_adapt==3 && mod(i,2) == 0)
          [xy,triQ,bks,ndone4]=adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options); 
          save_plot(tri,xy,(i-1)*4+4,options,sprintf('i=%0.0f, after smoothing',i));
       end;
    end
    if options.verbose == 2
  	     disp(sprintf('%2.0f%% coarsen, %2.0f%% flip, %2.0f%% refine, %2.0f%% move',ndone1*100,ndone2*100,ndone3*100,ndone4*100));
    end;
    if ismember(sw,[1,3,5,7,9,11,13,15])
    [xy,triQ,bks,ndone4]=adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options,5);
    save_plot(tri,xy,4*i+1,options,'after final smoothing');
    end
    if options.verbose
    	if nargin==9; disp('Recalculation qualities as angles in degrees');	options.qualP=eps;
    		triQ=elem_angle(tri,xy,options);
    	end;
    	disp_qual(triQ,tri,xy);
    end;
    if options.qualP < 0 && nargin == 8
    	options.qualP=-options.qualP; options.spltc=1; options.smpRFN=1; options.consRFN=2;
    	options.fastRFN=0;	options.consRM=1; options.fastRM=0;
    	disp(sprintf('Fixing angles larger than %0.1f degress', options.qualP*180/pi));
        triQ=elem_angle(tri,xy,options); disp_qual(triQ,tri,xy);
   	   	for i=1:options.MVit
   		[xy,tri,bndmesh,Nmetric,triQ,bks,ndone3]=adapt_add_nd(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
   		disp_qual(triQ,tri,xy,i)
        end;
    end;
  end


function [tri,xy,bndmesh,Nmetric,triQ,bks] =newadapt(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options)
    options.consRM=true;
    for i=1:options.innerit
               if options.verbose
	       disp_qual(triQ,tri,xy,i);
	end;
        while true
        	   if options.RMedg
        	   	[xy,tri,bndmesh,Nmetric,triQ,bks,ndone1]=adapt_rm_edg(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
        	   else
	            [tri,triQ,bks,bndmesh,ndone1,xy,Nmetric]=adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options);
	   end;
            if ndone1 < 1e-2
                break
            end;
        end;
        [xy,tri,bndmesh,Nmetric,triQ,ndone3]=adapt_add_nd(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
        while true %for j=1:3
            [xy,triQ,bks,ndone4]=adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options); 
            if not(size(xy,2) == 3 && options.swap3D == 0)
            [tri,triQ,bks,ndone2]=adapt_flipedg(tri,xy,Nmetric,bndmesh,triQ,bks,options); 
            else
            ndone2=0;
            end;
            if  size(xy,2) == 3
	    [tri,triQ,bks,bndmesh,ndone2_]=adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options);
	    ndone2=ndone2+ndone2_;
 	    end;
            if ndone2 < 1e-2 && ndone4 < 1e-2
                break
            end;
        end;
          if options.verbose
    	disp_qual(triQ,tri,xy);
     end;
        if ndone3  == 0 %&& ndone1 == 0 && ndone2 == 0
            disp('Adaptation has converged'); break;
        end;
    end;

function disp_qual(triQ,tri,xy,ii)
if nargin == 3
disp(sprintf('%6.0f nodes, %5.0f elements, %0.1f%% average quality (min=%2.1f%%))',size(xy,1),size(tri,1),mean(triQ)*100,min(triQ)*100));
else
disp(sprintf('%3.0f: %6.0f nodes, %5.0f elements, %0.1f%% average quality (min=%2.1f%%))',ii,size(xy,1),size(tri,1),mean(triQ)*100,min(triQ)*100));
end;

function save_plot(tri,xy,i4,options,intit)
if isfield(options,'forpres')
	if not(exist(options.forpres,'dir'))
		mkdir(options.forpres);
	end;
	if size(xy,2) == 2
	i4=i4-1;
	figure(); trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','none'); xlim([-0.01 1.01]); ylim([-0.01 1.01]); axis('off'); axis('equal'); box('off');  title(intit,'fontsize',32);
	print([options.forpres sprintf('/%0.0f.png',i4)],'-dpng');
	close all;
	else %3D
        if numel(strfind(version,'R2010a'))~=0
            [fac,fac2tri,tri2fac]=bks_all3D(tri);
            trimesh(fac(fac2tri(:,2)==0,:),xy(:,1),xy(:,2),xy(:,3),'edgecolor','k','marker','none'); axis off; box off;
            title(intit,'fontsize',32); print('-dpng',sprintf('./%s/%0.0f.png',options.forpres,i4));
        else
            save([options.forpres sprintf('/%0.0f.mat',i4)],'tri','xy','intit');
        end;
	end;
end;