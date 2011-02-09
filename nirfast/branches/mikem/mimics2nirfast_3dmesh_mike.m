function [mesh] = mimics2nirfast_3dmesh_mike(fn_mimics_mesh, fn_nirfast_mesh)

% [mesh] = mimics2nirfast_3dmesh_mike(fn_mimics_mesh, fn_nirfast_mesh)
%%Created Subha 2/10/09
%%This function creates element, node files in nirfast format from mimics
%%volume mesh exported in ansys format
%%Input filename is to be given without extension

% Modified s davis 2/11/09
%       Added param file
%       Automated loading of material properties by searching file for "!
%       Volume 1" and loading from that point on.


disp('Extracting nodes');
temp = importdata([fn_mimics_mesh,'.cdb'],' ');
nodes = temp.data;
nodes(:,1) = 0;


disp('Extracting elements');
clear temp;
fid = fopen([fn_mimics_mesh,'.cdb'],'r');
ln = 'junk'; it = 1;
while strcmp(ln,'(19i8)')~=1
    ln = fgetl(fid);
    it = it+1;
end
fclose(fid);
temp = importdata([fn_mimics_mesh,'.cdb'],' ',it-1);
a = find(temp.data(:,3)==1);
elem = temp.data(a,end-3:end);  
nel = length(elem);
clear it ln a
%reg = temp.data(a,1);
%reg(end) = temp.data(nel+1,1);


disp('Extracting surface nodes');
clear temp;
fid = fopen([fn_mimics_mesh,'.cdb'],'r');
ln = 'junk'; it = 1;
while strcmp(ln,'(1i10)')~=1
    ln = fgetl(fid);
    it = it+1;
end
fclose(fid);
surf_temp = importdata([fn_mimics_mesh,'.cdb'],' ',it-1);
surf_nodes = surf_temp.data;
nodes(surf_nodes,1) = 1; clear it ln


disp('Writing .node and .elem files');
fid = fopen([fn_nirfast_mesh,'.node'],'w');
for ii = 1:length(nodes)
    fprintf(fid,'%d %f %f %f \n',nodes(ii,1),nodes(ii,2),nodes(ii,3),nodes(ii,4));
end
fclose(fid);

fid = fopen([fn_nirfast_mesh,'.elem'],'w');
for ii = 1:nel
    fprintf(fid,'%d %d %d %d \n',elem(ii,1),elem(ii,2),elem(ii,3),elem(ii,4));
end
fclose(fid);

% Create region file
if exist([fn_mimics_mesh,'.txt'])~=0
    disp('Writing .region file from Mimics Materials list');
    
    % Check file and find Volume list header
    fid = fopen([fn_mimics_mesh,'.txt'],'r');
    ln = 'junk'; it = 1;
    while strcmp(ln,'! Volume 1')~=1
        ln = fgetl(fid);
        it = it+1;
    end
    fclose(fid);
    
    % Import data from Volume header on.  Extract material number.
    mat = importdata([fn_mimics_mesh,'.txt'],',',it);
    mat = mat.data;
    elem = elem(1:length(mat),:);
    region = zeros(length(nodes),1);
    for ii = 1:nel-1
        for jj = 1:4
            ind = elem(ii,jj);
            region(ind) = mat(ii,1);
        end
    end
    
    % Write region file
    fid = fopen([fn_nirfast_mesh,'.region'],'w');
    for ii = 1:length(region)
        fprintf(fid,'%d \n',region(ii));
    end
    fclose(fid);
    mesh.region = region; clear region
end

% Create param file
str = input('Please enter mesh type (1 = stnd, 2 = fluor, 3 = spec):');
fid = fopen([fn_nirfast_mesh,'.param'],'w');
if str == 1
    disp('Creating Param file for stnd mesh')
    fprintf(fid,'stnd\n');
    for ii = 1:length(nodes)
        fprintf(fid,'%f %f %f \n',0.01,0.330033,1.33);
    end
elseif str == 2
    disp('Creating Param file for fluor mesh')
    fprintf(fid,'fluor\n');
    for ii = 1:length(nodes)
        fprintf(fid,'%f %f %f %f %f %f %f %f \n',0.01,0.330033,1.33,0.01,0.330033,0.001,0.05,0);
    end
elseif str == 3
    disp('Creating Param file for spec mesh')
    fprintf(fid,'spec\n');
    fprintf(fid,'HbO\n');
    fprintf(fid,'deoxyHb\n');
    fprintf(fid,'Water\n');
    fprintf(fid,'S-Amplitude\n');
    fprintf(fid,'S-Power\n');
    for ii = 1:length(nodes)
        fprintf(fid,'%f %f %f %f %f \n',0.01,0.01,0.5,1,1);
    end
end
fclose(fid);

mesh.nodes = nodes; clear nodes
mesh.elements = elem; clear elem