function mesh = move_detector(mesh)

% mesh = move_detector(mesh)
%
% Moves the detectors onto the surface of the mesh
%
% mesh is the mesh location or variable


%% load mesh
if ischar(mesh)
  mesh = load_mesh(mesh);
end
if ~isfield(mesh,'meas') || ~isfield(mesh.meas,'coord')
    errordlg('No detectors present','NIRFAST Error');
    error('No detectors present');
end

if size(mesh.meas.coord,2) == 2
    mesh.meas.coord(:,end+1) = 0;
end

%% get list of boundary faces
if size(mesh.elements,2) == 4
    faces = [mesh.elements(:,[1,2,3]);
              mesh.elements(:,[1,2,4]);
              mesh.elements(:,[1,3,4]);
              mesh.elements(:,[2,3,4])];
    faces = sort(faces,2);
    faces = unique(faces,'rows');
    faces = faces(sum(mesh.bndvtx(faces),2)==3,:);
elseif size(mesh.elements,2) == 3
    if mesh.dimension == 3
        faces = mesh.elements(sum(mesh.bndvtx(mesh.elements),2)==3,:);
    elseif mesh.dimension == 2
        faces = mesh.elements(sum(mesh.bndvtx(mesh.elements),2)==2,:);
    end
end

%% loop through detectors
for i=1:size(mesh.meas.coord,1)
    
    % find closest boundary node
    if mesh.dimension == 2
        dist = distance(mesh.nodes,mesh.bndvtx,[mesh.meas.coord(i,:) 0]);
    elseif mesh.dimension == 3
        dist = distance(mesh.nodes,mesh.bndvtx,mesh.meas.coord(i,:));
    end
    r0_ind = find(dist==min(dist));
    r0_ind = r0_ind(1);
    
    % find faces including the closest boundary node
    fi = faces(sum(faces==r0_ind,2)>0,:);

    % find closest face
    dist = zeros(size(fi,1),1);
    point = zeros(size(fi,1),3);
    for j=1:size(fi,1)
        [dist(j),point(j,:)] = pointTriangleDistance([mesh.nodes(fi(j,1),:);...
            mesh.nodes(fi(j,2),:);mesh.nodes(fi(j,3),:)],mesh.meas.coord(i,:));
    end
    smallest = find(dist == min(dist));

    % move detector to the closest point on that face
    if size(mesh.meas.coord,2) == 2
        mesh.meas.coord(i,:) = point(smallest(1),1:2);
    elseif size(mesh.meas.coord,2) == 3
        mesh.meas.coord(i,:) = point(smallest(1),:);
    end
        
end