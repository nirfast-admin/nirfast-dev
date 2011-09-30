function nirfast2vtk(meshfname,listsolfnames,iter)
% This function takes a NIRFAST mesh and an associated list of solution 
% files to produce a combined .vtk file for viewing/manipulation in
% Paraview and other visualization software. Filenames should be provided 
% without extensions. 
%
% usage: nirfast2vtk('mesh_fn',{fname_HbO,fname_Water...},8);
% outputs: 'mesh_fn_wsol_iter8.vtk'
%
% author: venkat krishnawamy/03242010
% last update: 
% part of NIRFAST package
% (C) Hamid Deghani 2008

outfname = [meshfname,'_wsol_iter',num2str(iter),'.vtk'];
numsol = length(listsolfnames);

% read-in representation from NIRFAST files
temp = dlmread([meshfname '.node']);
nodes = temp(:,end-2:end);
numnodes = length(nodes);
elems = dlmread([meshfname '.elem']);
numelems = length(elems);

%read-in solutions from NIRFAST .sol files
soldata = zeros(numnodes,numsol);
for i = 1:numsol
    fid = fopen([listsolfnames{i}, '.sol']);
    temp = textscan(fid,'%f ','HeaderLines', (iter-1)*2+1);
    soldata(:,i) = temp{:};
    fclose(fid);
end;

fid = fopen(outfname,'w');

%define VTK header for FEM mesh representation
line0 = '# vtk DataFile Version 2.0';
line1 = ['NIRFAST mesh with solutions for case: {', meshfname,'}'];
line2 = 'ASCII';
line3 = 'DATASET UNSTRUCTURED_GRID';
fprintf(fid,'%s\n%s\n%s\n',line0,line1,line2,line3);

line4 = ['POINTS ', num2str(numnodes), ' float']; %node defs
fprintf(fid,'%s\n',line4);
fprintf(fid, '%f %f %f\n', nodes');

line5 = ['CELLS ',num2str(numelems),' ',num2str(numelems*4+numelems)]; %connectivity maps
fprintf(fid,'%s\n',line5);
fprintf(fid,'%d %d %d %d %d\n',[4*ones(numelems,1) elems-1]');
line6 = ['CELL_TYPES ', num2str(numelems)]; %specify the mesh basis 10-tetrahedral for all connectivity maps 
fprintf(fid,'%s\n',line6);
fprintf(fid,'%d\n', ones(numelems,1)*10);
fprintf(fid,'%s\n',['POINT_DATA ',num2str(numnodes)]); %specify the data that follows is defined on the nodes

for i = 1:numsol
    fprintf(fid,'%s\n',['SCALARS ', listsolfnames{i}, ' float 1']);
    fprintf(fid,'%s\n','LOOKUP_TABLE default');
    fprintf(fid,'%f\n', soldata(:,i));
end;

fclose(fid);