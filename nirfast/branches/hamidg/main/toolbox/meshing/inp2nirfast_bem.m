function inp2nirfast_bem(fnprefix,saveloc,source_coord_fn)

% inp2nirfast_bem(fnprefix,saveloc,source_coord_fn)
%
% Converts inp files to a nirfast bem mesh
%
% fnprefix is the prefix location of the inp files
% saveloc is the location to save the mesh to
% source_coord_fn is a text file containing coordinates of the source
% locations. If not provided, this routine will atempt to find
% 'fnprefix_measurements_1.txt' and use the data from that. 


no_regions = length(dir([fnprefix '*.inp']));
if no_regions==0
    error(['Cannot find file .inp files whose prefix is ' fnprefix]);
end

fprintf('\n\tConverting inp files and re-orienting\n');
for i=1:no_regions
    fn = [fnprefix num2str(i) '.inp'];
    abaqus2nodele_surface(fn);
end
fprintf('\tSurface detection.\n\t\t');
surface_relations_mex(fnprefix,no_regions);
relations = textread('surface_relations.txt');

region_id=zeros(no_regions,1);
idcounter=1;
for i=1:size(relations,1)
    if i==1
        fn = [fnprefix num2str(relations(i,1))];
        [ele,node] = read_nod_elm(fn,1);
        region1 = repmat(1,[size(ele,1) 1]);
        region2 = repmat(0,[size(ele,1) 1]);
        no_ext_nodes = size(node,1);
        NoBdyNodes(i) = no_ext_nodes;
        NoBdyElems(i) = size(ele,1);
    end
    if region_id(relations(i,1)) == 0
        region_id(relations(i,1)) = idcounter;
        idcounter = idcounter + 1;
    end
    momid = region_id(relations(i,1));
    for j=2:size(relations,2)
        if relations(i,j)==0, continue; end
            fn = [fnprefix num2str(relations(i,j))];
            kidid = idcounter;
            region_id(relations(i,j)) = idcounter;
            
            [fooele,foonode] = read_nod_elm(fn,1);
            ele=[ele;fooele+size(node,1)]; %#ok<AGROW>
            node=[node; foonode]; %#ok<AGROW>
            region1 = cat(1,region1,repmat(momid,[size(fooele,1) 1]));
            region2 = cat(1,region2,repmat(kidid,[size(fooele,1) 1]));
            
            NoBdyNodes(idcounter) = size(foonode,1);
            NoBdyElems(idcounter) = size(fooele,1);
            idcounter = idcounter + 1;
    end
end
foo=zeros(size(node,1),1);
foo(1:no_ext_nodes,:) = 1;
node=[foo node];

fprintf('\tWriting nirfast mesh files\n');
dlmwrite([saveloc '.node'],node,' ');
dlmwrite([saveloc '.elem'],ele,' ');
dlmwrite([saveloc '.region'],[region1 region2],' ');

if nargin < 3
% try to find 'fnprefix_measurements_1.txt' that can contain source
% locations exported by mimics
    if fnprefix(end)=='_', fnprefix(end)=[]; end
    foo = dir([fnprefix '_measurements_*.txt']);
    if length(foo)==1
        source_coord_fn = foo.name;
    else
        source_coord_fn = [];
    end
end
if ~isempty(source_coord_fn)
    scoords = importdata(source_coord_fn);
    if isstruct(scoords)
        scoords=scoords.data(:,2:4);
    end
    dlmwrite([saveloc '.source'],scoords,'delimiter',' ','precision','%0.12f');
end




