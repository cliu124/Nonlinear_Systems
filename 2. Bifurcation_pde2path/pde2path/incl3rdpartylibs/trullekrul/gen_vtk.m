function [] = gen_vtk(flname,tri,xy,sclrvr,ii)
%tabchar = '\t'

fid = fopen(flname,'w');
output1 = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n  <Collection>\n';
output2 = [];
for i=1:ii
output2 = [output2 sprintf('    <DataSet timestep="%d" part="0" file="%s%d.vtu" />\n', i-1, flname(1:end-4), i-1)];
%output2 = sprintf('<DataSet timestep="%d" part="0" file=".vtu" />\n'., i);
end;% output2 = output2(1:end-2);
output3 = '  </Collection>\n</VTKFile>';
output = [output1 output2 output3];
fprintf(fid,output);
fclose(fid);

tabchar = '  ';
fid = fopen(sprintf('%s%d.vtu',flname(1:end-4),ii-1),'w');
output1 = sprintf('<?xml version="1.0"?>\\n<VTKFile type="UnstructuredGrid"  version="0.1"  >\\n<UnstructuredGrid>\\n<Piece  NumberOfPoints="%d" NumberOfCells="%d">\\n<Points>\\n<DataArray  type="Float64"  NumberOfComponents="%d"  format="ascii">',size(xy,1),size(tri,1),3); %size(xy,2)
if size(xy,2) == 2
 output2 = sprintf(['%6e %6e %6e' tabchar],[xy zeros(size(xy,1),1)]');
 output4 = sprintf(['%d %d %d' tabchar],tri'-1); %-repmat((1:size(tri,1))',1,3)
else %==3
 output2 = sprintf(['%6e %6e %6e' tabchar],xy');
 output4 = sprintf(['%d %d %d %d' tabchar],tri'-1);
end;
%output6 = sprintf('%d ',zeros(size(tri,1),1));
output6 = sprintf('%d ',size(tri,2):size(tri,2):numel(tri));
output8 = sprintf('%d ',(5+5*(size(xy,2)==3))*ones(size(tri,1),1));
output10 = sprintf(['%6e' tabchar], sclrvr);
output3 = '</DataArray>\n</Points>\n<Cells>\n<DataArray  type="UInt32"  Name="connectivity"  format="ascii">';
output9 = '</DataArray>\n</Cells>\n<PointData  Scalars="scalar variable">\n<DataArray  type="Float64"  Name="scalar variable"  format="ascii">';
output11 = '</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>';
output5 = '</DataArray>\n<DataArray  type="UInt32"  Name="offsets"  format="ascii">';
output7 = '</DataArray>\n<DataArray  type="UInt8"  Name="types"  format="ascii">';
output = [output1 output2 output3 output4 output5 output6 output7 output8 output9 output10 output11];
fprintf(fid,output);
fclose(fid);