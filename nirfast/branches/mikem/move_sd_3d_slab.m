function mesh = move_sd_3d_slab(mesh,w,fn_oppfibers)

% This function takes the s-d mesh *.source and *.meas files created using
% mimics2nirfast_sd_3D_mouse and moves the source and detector locations to
% be compatible with Nirfast modeling.  Specifically, sources are moved 1
% scattering distance in and detectors are moved just inside the mesh
% surface.
% CAUTION:  This program is for radially arrays of fibers,
% commonly used in MRI mouse imaging in 3-D.  It will not accomodate slab
% geometries or 2-D information.
% This program assumes sources and detectors are colocalized and # sources
% = # detectors

% Inputs:
% mesh is the mesh filename or workspace variable complete with mesh.source
% and mesh.meas information produced from mimics2nirfast_sd_3D_mouse
% w is width of gaussian source


if ischar(mesh)==1
    mesh = load_mesh(mesh);
end

% Make sure first coolumn of opposing fibers list is 1:#sources, matching
% the source coord.
opp_fibers = load(fn_oppfibers);

%****************************************************
% Move s-d positions to surface of mesh:

% move opposing fibers toward one another and place on mesh surface:
[nsource,junk] = size(mesh.source.coord);
% error check
if length(opp_fibers)~=nsource
    disp('Error: length of opposing fibers list does not match length of source.coord')
    return
end

for i = 1 : nsource
    % distance of each boundary node from source
    dist = distance(mesh.nodes,mesh.bndvtx,[mesh.source.coord(i,:) 0]);
    
    % sort nodes in order of nearest to farthest
    [snodes,ind] = sort(dist);
    
    % Start with nearest node and test surface elements of which the node is a
    % vertex
    j=1; true = 0;
    
    % find surface elements with node 'x' as a vertex
    while true == 0
        % ind is list of node numbers in order of distance from
        % source/det point.  Find the elements which contain this node:
        [r,c,v] = find(mesh.elements == ind(j));
        
        % foo is matrix of all nodes connected to node 'x', in node
        % connectivity list form.  We'll call it a 'shortened connectivity list'.
        % Contains both surface and internal nodes.
        foo = mesh.elements(r,:);
        [n,m] = size(foo);
        foo = reshape(foo,numel(foo),1);
        
        % set non-boundary nodes to NaN in shortened connectivity list:
        for k = 1:numel(foo)
            if mesh.bndvtx(foo(k))==0;
                foo(k) = NaN;
            end
        end
        foo = reshape(foo, n, m);
        
        % For any row in shortened connectivity list, check to see if
        % the number of NaN = 1.  These rows represent the boundary
        % elements
        foo2 = [];
        for k = 1:n
            nans = find(isnan(foo(k,:)));
            if numel(nans) == 1
                temp = foo(k,:);
                temp(find(isnan(temp)==1))=[];
                foo2 = [foo2 ; temp];  %foo2 is matrix of surface node numbers connected to node 'x'
            end
        end
        clear foo
        
        % Test surface elements to see if line between opposing s-d
        % pairs intersects with surface element.  If it intersects a particular element,
        % set the intersection point as the new s-d coordinate.
        [n,m] = size(foo2);
        k = 1; true = 0;
        while true == 0 & k <= n
            % To make the syntax a little easier to read, define points P, Q, R - vertices of the surface triangle
            % which we are testing
            P = mesh.nodes(foo2(k,1),:); Q = mesh.nodes(foo2(k,2),:); R = mesh.nodes(foo2(k,3),:);
            
            % fit plane defined by P, Q, R
            plane = fitplane([P; Q; R]');
            
            % determine intersection of line segment between opposing
            % s-d pair with plane PQR
            [I,check]=plane_line_intersect(plane(1:3),P,mesh.source.coord(i,:),mesh.source.coord(opp_fibers(i,2),:));
            
            % check to see if intersection point is within element:
            true=inside_triangle(I,P,Q,R);
            if true == 1
                mesh.source.coord(i,:) = I;
            elseif true == 0
                k = k+1;
            end
        end
        j=j+1;
    end
    % Plot points to see if fiber location is inside surface element.
    %     figure
    %     xx = [P(1);Q(1);R(1)];
    %     yy = [P(2);Q(2);R(2)];
    %     zz = [P(3);Q(3);R(3)];
    %     plot3(xx,yy,zz,'-O',I(1),I(2),I(3),'rx')
end
mesh.meas.coord = mesh.source.coord;

%***********************************************
% Now move sources in 1 scatter distance from surface and detectors just
% inside mesh (by 0.001)
if strcmp(mesh.type,'fluor')==1
    mus_eff = mesh.musx;
elseif strcmp(mesh.type,'stnd')==1
    mus_eff = mesh.mus;
elseif strcmp(mesh.type,'spec')==1
    [mua, mus_eff, kappa, E] = calc_mua_mus(mesh,785);
    clear mua kappa E
end

for i = 1:nsource
    % distance of each boundary node from source
    dist = distance(mesh.nodes,mesh.bndvtx,mesh.source.coord(i,:));
    % index of nearest boundary node
    r0_ind = find(dist==min(dist));
    r0_ind = r0_ind(1);
    % mean scatter value of where source will be
    dist = distance(mesh.nodes,ones(length(mesh.bndvtx),1),...
        [mesh.nodes(r0_ind,1) ...
        mesh.nodes(r0_ind,2) ...
        mesh.nodes(r0_ind,3)]);
    scat_dist = 1./mean(mus_eff(find(dist<=w)));
    
    fiber_vec = (-mesh.source.coord(i,:)+mesh.source.coord(opp_fibers(i,2),:))...
        /norm(-mesh.source.coord(i,:)+mesh.source.coord(opp_fibers(i,2),:));
    
    % move sources
    mesh.source.coord(i,:) = [fiber_vec*scat_dist+mesh.source.coord(i,:)];
    
    % move detectors, check int functions and move further if needed
    movein = 0.001;
    mesh.meas.coord(i,:) = [fiber_vec*movein+mesh.meas.coord(i,:)];
%     [junk,int_func(i,:)] = tsearchn(mesh.nodes,...
%         mesh.elements,...
%         mesh.meas.coord(i,:));
%     while isnan(int_func(end,:)) == 1
%         movein = 10*movein;
%         mesh.meas.coord(i,:) = [fiber_vec*movein+mesh.meas.coord(i,:)];
%         [junk,int_func(i,:)] = tsearchn(mesh.nodes,...
%             mesh.elements,...
%             mesh.meas.coord(i,:));
%    end
end
[ind,int_func] = tsearchn(mesh.nodes,...
    mesh.elements,...
    mesh.meas.coord);
mesh.meas.int_func = [ind int_func]


%*******************************************
% Finally, fix calculated positions so s-d pairs are not moved when loaded
% using load_mesh
mesh.meas.fixed = 1;
mesh.source.fixed = 1;
