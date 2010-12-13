function [data,u,mesh]= bemdata_stnd(mesh,frequency,myargs)
% [data,u,mesh]= bemdata_stnd(mesh,frequency)
%
% Calculates data (phase and amplitude) for a given
% mesh at a given frequency (MHz).
% outputs phase and amplitude in structure data
% and mesh information in mesh
% myargs: used to pass additional information such --verbose flag
% 'data.ppa' is the phase and amplitude values at source/detector locations
% 'data.phi' is the field value for all nodes of the boundaries.
% Written By:
%           Hamid R Ghadyani, March 2010

%% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

%% load mesh
if ischar(mesh)
    mesh = load_mesh(mesh);
end
if isfield(mesh,'region')==0
    errordlg([mesh.name ' mesh needs to have a .region field!'],'NIRFAST Error');
    error([mesh.name ' mesh needs to have a .region field!']);        
end

verbose=0;
if nargin==3
    if isfield(myargs,'verbose')
        verbose=myargs.verbose;
    end
end

c=(3e11./mesh.ri);
omega = sqrt((mesh.mua + (1i*(2*pi*frequency*1e6)./c))./mesh.kappa);
% Establish boundary relations based on .region file.
relations = GetSurfRelations(mesh);
mesh. relations = relations;
[NoBdyNodes regionInfo BdyNodeIndices]=GetNoBdyNodes(mesh);
% Initialze RHS index vector
totn = NoBdyNodes(1);
for i=2:length(NoBdyNodes)
    totn = totn + 2*NoBdyNodes(i);
end
rhs_idx = false(totn,1);
q_tot=zeros(totn,1);
K=zeros(totn,totn);

num_sources = size(mesh.source.coord,1);
nnod = length(mesh.nodes);
visits=ones(size(relations,1),1);
%% Main loop (over all the sources)
for scounter=1:num_sources
    tic
    if scounter==1
        if verbose==1, fprintf('.'); elseif verbose==2, fprintf('    Building components of BEM matrix for region %d... ',rid); end
        % Construct the LHS matrix 'K' only once
        for region=1:size(mesh.relations,1);
            rid=relations(region,1);
            [region_elems region_nodes] = GetNodesAndElementsOfRegion(mesh,regionInfo(rid));
            %region_coords = mesh.nodes(region_nodes,:);
            
            [ar ai br bi nids] = BuildMainMatrix_K(mesh,region_elems,omega(region),mesh.kappa(region),16);
            %[ar ai br bi] = main_build_matrix_K(mesh.nodes, region_elems, region_coords, region_nodes,...
                %omega(region), mesh.kappa(region), 2048);
            A = complex(ar,ai); B = complex(br,bi);
            
            % Take care of interior solid angle terms
            tmp=length(A);
            idx=1:tmp+1:tmp^2;
            %A(idx(region_nodes)') = 1 - (sum(A(region_nodes,region_nodes),2) - A(idx(region_nodes)'));
            A(idx') = 1 - (sum(A,2) - A(idx'));

            foor = relations(region,:); foor = setdiff(foor,0);
            for i=1:length(foor)
                Node_idxm1 = BdyNodeIndices{foor(i)};
                [tf Node_idxm] = ismember(Node_idxm1,nids);
                for j=1:length(foor)
                    Node_idxn1 = BdyNodeIndices{foor(j)};
                    [tf Node_idxn] = ismember(Node_idxn1,nids);
                    % Get locations of current sub-matrix in big K
                    [mstart mend nstart nend] = CalculateLocationInK(...
                        foor(i),foor(j),NoBdyNodes,visits(foor(i)));
                    % Save the location of nodes of region I (one) to apply
                    % BC's
                    if rid==1 && j==1
                        rhs_idx(mstart:mend) = true(mend-mstart+1,1);
                    end
                    if foor(j)==1
                        K(mstart:mend,nstart:nend) = [A(Node_idxm,Node_idxn) + (mesh.ksi(rid)).*B(Node_idxm,Node_idxn)];
%                         K(mstart:mend,nstart:nend) = [A + (mesh.ksi(rid)).*B];
                    elseif foor(j)~=rid
                        K(mstart:mend,nstart:nend) = [A(Node_idxm,Node_idxn)  -B(Node_idxm,Node_idxn)];
%                         K(mstart:mend,nstart:nend) = [A  -B];
                    else
                        K(mstart:mend,nstart:nend) = [A(Node_idxm,Node_idxn) B(Node_idxm,Node_idxn)];
%                         K(mstart:mend,nstart:nend) = [A B];
                    end
                end
            end
            bf=relations(region,:)~=0;
            visits(relations(region,bf)) = visits(relations(region,bf)) + 1;
        end
        if verbose==1, fprintf('.\n'); elseif verbose==2, fprintf('done!\n'); end
        clear ar ai br bi A B
        % Store nodes and elements of region I
        [regionI_elems regionI_nodes] = GetNodesAndElementsOfRegion(mesh,regionInfo(relations(1,1))); 
    end
    toc
    qq = build_source(nnod,mesh.nodes(regionI_nodes,:),regionI_nodes,...
        omega(1),mesh.kappa(1),mesh.source.coord(scounter,:),1);
    q_tot(rhs_idx,scounter) = qq(regionI_nodes,:);
end
% Solve the equation: K.u = q_tot
tic
u = K\q_tot;
toc
% Get Nodal solutions and phase/amplitude at detector locations
u = GetNodalSolutions(u,mesh,NoBdyNodes);
data=GetDataAtDetectorLocations(u,mesh);





















function [ar ai br bi nids] = BuildMainMatrix_K(mesh,region_elems,omega,D,num_procs)
% Call the MPI-based matrix build and wait till it's done.
warning('off','MATLAB:DELETE:FileNotFound');

% Renumber the region nodes and elements such that nodes are from 1 to
% number of nodes in region
region_elems=region_elems(:,2:4);
nids = unique([region_elems(:,1);region_elems(:,2);region_elems(:,3)]);
[tf elements]=ismember(region_elems,nids);
nodes = mesh.nodes(nids,:);

save('bem_mpi_input.mat','nodes','elements','omega','D','num_procs');
delete('bem_mpi.lock');
command = which('bem_build_main_matrix_mpi.exe');
mpidir = fileparts(command);
syscommand=['mpirun -np ' num2str(num_procs) ' -machinefile ' mpidir '/machinefile ' command];
system(syscommand);
while true
    if exist('bem_mpi.lock','file')
        break
    end
end
warning('on','MATLAB:DELETE:FileNotFound');
delete('bem_mpi.lock');
result = load('Ar.mat'); ar = result.Ar;
result = load('Ai.mat'); ai = result.Ai;
result = load('Br.mat'); br = result.Br;
result = load('Bi.mat'); bi = result.Bi;
delete('Ar.mat','Ai.mat','Br.mat','Bi.mat');




function data = GetDataAtDetectorLocations(u,mesh)
% Calculates phase and amplitude values at detector locations
%%
tmp=[];
for scounter = 1:size(mesh.source.coord,1)
    detectors_used = mesh.link(scounter,:);
    detectors_used = detectors_used(find(detectors_used(:) ~= 0)');
    meas_int_func_ss = mesh.meas.int_func(detectors_used',:);
    nrr = size(meas_int_func_ss,1);
    if (nrr == 0)
        data_ss = [];
    else
        data_ss = get_boundary_data_bem(mesh.elements,meas_int_func_ss,u(:,scounter));
    end
    tmp = [tmp; data_ss'];
end

data.paa = [abs(tmp) angle(tmp).*180/pi];
data.phi = u;

function u = GetNodalSolutions(u,mesh,NoBdyNodes)
% Drop solutions for flux values and just return the nodal intensity solutions
%%

relations=mesh.relations;
startidx = NoBdyNodes(1);
bf=false(size(u,1),1);
bf(1:startidx) = true;

endidx=startidx;
for i=2:size(relations,1)
    startidx = endidx + 1;
    endidx = startidx+NoBdyNodes(i)-1;
    bf(startidx:endidx)=true;
    endidx = endidx + NoBdyNodes(i);
end
assert(sum(bf)==sum(NoBdyNodes));

u=u(bf,:);

function [mstart mend nstart nend] = CalculateLocationInK(m,n,NoBdyNodes,visitn)
% Calculate the location of sub-matrices that build our main A and B
% matrices.
%%

if m==1
    mstart = 1;
    mend = NoBdyNodes(m);
else
    mstart = NoBdyNodes(1);
    foo = 0;
    for i=2:(m-1)
        foo = foo + 2*NoBdyNodes(i);
    end
    mstart = mstart + foo + 1 + (visitn-1)*NoBdyNodes(m);
    mend = mstart + NoBdyNodes(m) - 1;
end

if n==1
    nstart = 1;
    nend = NoBdyNodes(n);
else
    nstart = NoBdyNodes(1);
    foo = 0;
    for i=2:(n-1)
        foo = foo + 2*NoBdyNodes(i);
    end
    nstart = nstart + foo + 1 ;
    nend = nstart + 2*NoBdyNodes(n) -1;
end

function [NoBdyNodes regionInfo BdyNodeIndices]=GetNoBdyNodes(mesh)
% Calculate how many nodes each boundary has
% We assume that the most exterior shell has boundary ID of 1 and its
% mesh.region looks like (1,0). The other boundaries will have mesh.region 
% as (id1,id2) where id1 is the boundary that encloses boundary with id2
%%

relations=mesh.relations;
NoBdyNodes=zeros(size(relations,1),1);
BdyNodeIndices=cell(size(relations,1),1);

for i=1:size(relations,1)
    region=relations(i,1);
    % find all the region IDs involved with 'region'
    allflag=mesh.region(:,1)==region | mesh.region(:,2)==region;
    allregions = unique([mesh.region(allflag,1); mesh.region(allflag,2)]);
    allregions = setdiff(allregions,region);

    % find region IDs of direct children of 'region'
    myregions=[];
    for j=2:size(relations,2)
        if relations(i,j)==0, continue; end
        myregions=[myregions relations(i,j)];
    end

    % find the region ID that 'region' resides in
    myext_region = setdiff(allregions,myregions);
    if length(myext_region)~=1
        error('Corrupt Mesh! Check your mesh.region structure.');
    end
    regionInfo(region).allflag=allflag;
    regionInfo(region).myext_region=myext_region;
    regionInfo(region).region=region;
    
    for j=1:size(relations,2)
        if relations(i,j)==0, continue; end
        [Node_idxm] = GetNodesForBdy(mesh,relations(i,j),regionInfo(region));
        NoBdyNodes(relations(i,j))=length(Node_idxm);
        BdyNodeIndices{relations(i,j)}=Node_idxm;
    end
end



function [relations]= GetSurfRelations(mesh)
% Using mesh.region, build relations matrix, which is basically telling us
% every region's immediate sub-regions.
%%

ids = unique([mesh.region(:,1); mesh.region(:,2)]);
ids=ids(2:end); % Removing '0'

% Build 'surface_relations' matrix
for i=1:length(ids)
    bf = mesh.region(:,1)==ids(i);
    inter_regions = unique(mesh.region(bf,2));
    [tf idx]=ismember(0,inter_regions);
    if tf
        inter_regions(idx) = [];
    end
    inter_regions = setdiff(inter_regions,ids(i));
    if isempty(inter_regions), inter_regions=[0]; end;
    foo=[ids(i) inter_regions'];
    relations(i,1:length(foo)) = foo;
end


function [Node_idxm] = GetNodesForBdy(mesh,m,regionInfo)
% Returns the indecies of all the nodes that belong to boundary 'm'
%%
myext_region = regionInfo.myext_region;
allflag = regionInfo.allflag;

region_elms = mesh.elements(allflag,:);
region_reg = mesh.region(allflag,:);

if m==regionInfo.region
    rid = myext_region;
else
    rid = m;
end
foo = region_elms(region_reg(:,1)==rid | region_reg(:,2)==rid,:);
Node_idxm = unique(reshape(foo,[],1));


function [region_elems region_nodes] = GetNodesAndElementsOfRegion(mesh,regionInfo)
% Returns elements and nodes that define 'region' in the 'mesh'
% All the elemnts that belong to the should have their
% normal vectors pointing outward. To do this, we assume
% that all the elements present in mesh.elements are
% already oriented and they are pointing outward with respect to their own
% interior.
% this function also returns 'myext_region' which is the region ID that our
% 'region' resides in

% get the region ID that 'region' resides in
myext_region = regionInfo.myext_region;
allflag = regionInfo.allflag;

% find the elements that are enclosing the rest
idx=(1:length(allflag))';idxmap=idx(allflag);
exterior_elms = mesh.elements(idxmap(mesh.region(allflag,1)==myext_region | mesh.region(allflag,2)==myext_region),:);

% find the interior elements and change their orientation so that their
% normal is pointing outward the current 'region'
other_elms = setdiff(mesh.elements(allflag,:),exterior_elms,'rows');
other_elms = [other_elms(:,1) other_elms(:,3) other_elms(:,2)];

% get all the nodes of the elements
region_elems = [exterior_elms; other_elms];
region_nodes = unique(reshape(region_elems,[],1));

region_elems = [(1:size(region_elems,1))' region_elems];
