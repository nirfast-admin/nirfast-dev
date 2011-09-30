function nirfast2tecplot(meshfname,listsolfnames)
% This function takes a NIRFAST mesh and associated solution files and
% produces a .dat file for viewing/manipulation in TecPlot. Input filenames
% should be provided without extensions.
%
% venkat krishnawamy/02102010
% part of tomo-tools package

outfname = [meshfname,'_sol_tecplot.dat'];
sollabels = '';

disp('reading mesh solution files...')

for i = 1:size(listsolfnames,2)
fid = fopen([listsolfnames{i}, '.sol']);
temp = textscan(fid,'%s %s \n');
sollabels = [sollabels, ', "', cell2mat(temp{2}),'"'];
temp = textscan(fid,'%f ','HeaderLines', 1);
soldata(:,i) = temp{:};
fclose(fid); 
end;

temp = dlmread([meshfname '.node']);
nodes = temp(:,end-2:end);
nodes_with_sol = [nodes soldata];
elems = dlmread([meshfname '.elem']);

line1 = ['TITLE = NIRFAST solution data file for case: {', meshfname,'}'];
line2 = ['VARIABLES = "X", "Y", "Z"' sollabels];
line3 = ['ZONE N = ',num2str(length(nodes)),' E = ',num2str(length(elems)),', DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON'];

disp('writing mesh and solutions data in Tecplot format...');

fid = fopen(outfname,'w');
fprintf(fid,'%s\n%s\n%s\n',line1,line2,line3);
fclose(fid);

dlmwrite(outfname, nodes_with_sol,'-append');
dlmwrite(outfname, elems,'-append');
disp('done!');