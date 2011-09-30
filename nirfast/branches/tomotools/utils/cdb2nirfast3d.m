function cdb2nirfast3d(fname, varargin)
% This converter is compatible only with ANSYS .cdb files generated using
% MIMICS v13.0. Slight modifications in the offsets are necessary to handle
% exports from other versions of MIMICS.
% 
% If containing multiple regions, the subvolume order must start from the
% largest e.g. whole (vol1), fibroglandular (vol2) and tumor (vol3)
%
% usage: cdb2nirfast3d(fname) for no region tagging
%        cdb2nirfast3d('TagRegions',true) for region tagging
%
% venkat krishnawamy/02012010
% part of tomo-tools package

p = inputParser;
p.addRequired('fname', @ischar);
p.addParamValue('TagRegions',false,@islogical);
p.parse(fname,varargin{:});

inpfname = [fname '.cdb'];
fid = fopen(inpfname);
temp = textscan(fid, '%f %f %f %f', 'HeaderLines', 7);
nodes = [temp{1} temp{2} temp{3} temp{4}];
nodes(:,1) = 0;

temp = textscan(fid,'%d %*f %*f %*f %*f %*f %*f %*f %*d %*f %*f %f %f %f %f', 'HeaderLines',7);
temp2 = [temp{1} temp{2} temp{3} temp{4} temp{5}];
temp2 = temp2(1:end-1,:); %get rid of that -1
numsurf = length(unique(temp2(:,1)));
elems = temp2(:, end-3:end);

regnodes = zeros(length(nodes),1);
% if a surface node is shared between two sub-volumes, one enclosing the
% other, tag it as belonging to the enclosed sub-volume.
if(p.Results.TagRegions)    
    for i = numsurf:-1:1
        regtemp = temp2(temp2(:,1) == i,end-3:end);
        regtemp = unique(regtemp);
        regnodes(regtemp,1)=i;
    end;
else    
end;

% only nodes lying on the outermost surface are tagged as 'surface nodes' in
% NIRFAST, inner boundaries are not handled for index mismatches
temp = textscan(fid, '%d', 'HeaderLines', 3);
outsurfnodes =  temp{1};
nodes(outsurfnodes,1) = 1;

params = repmat([0.01,0.33,1.33],length(nodes),1);

disp('saving nodes...')
dlmwrite([fname '.node'],nodes,'precision',6,'delimiter',' ');

disp('saving elems...')
dlmwrite([fname '.elem'],elems,'precision',6,'delimiter',' ');

disp('saving regions...')
dlmwrite([fname '.region'],regnodes,'precision',6,'delimiter',' ');

disp('saving standard parameters...')
dlmwrite([fname '.param'],params,'precision',6,'delimiter',' ');

disp('length of region data: '); disp(length(regnodes));
disp('length of node data: '); disp(length(nodes));

fclose(fid);

