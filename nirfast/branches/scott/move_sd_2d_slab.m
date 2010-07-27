function mesh = move_sd_2D_slab(mesh,w,fn_oppfibers)

% This function takes the s-d mesh *.source and *.meas files
% and moves the source and detector locations to
% be compatible with Nirfast modeling.  Specifically, sources are moved 1
% scattering distance in and detectors are moved just inside the mesh
% surface.  The assumption is that each fiber has an opposing fiber on the
% opposite side of the tissue, as specified by fn_oppfibers.  For examples,
% in a slab geometry, fibers on either side of the tissue volume face one
% another.  This programs moves them toward or away from their opposing
% fiber to the correct coordinate.

% CAUTION:  This program will not accomodate 2-D information.
% This program assumes sources and detectors are colocalized and # sources
% = # detectors.  Once the coordinates have been found, it is a simple
% matter to delete positions that are not actually being used (for example,
% some slab geometries are arranged with sources on one side of the tissue
% and dectectors on the other).

% Details:
% 1.  Loops over source fiber positions input in mesh.source.coord.
% 2.  For each source coord, orders boundary nodes based on distance to
% the coordinate (shortest to longest).
% 3.  Finds all neighboring surface nodes to the node under consideration.
% 4.  Each group of three surface nodes, including the node under
% consideration, forms a triangle on the surface of the mesh.  Loop over
% each triangle and do the following:
%       a.  Make a plane with the 3 points and a line segment with the
%       source coordinate and the coordinate of the opposing fibers
%       b.  Check to see if the line segment intersects the plane.  If so,
%       find the intersection point - this may well be outside of the
%       triangle of 3 points that define the plane.
%       c.  If there is an intersection point, check to see if it is inside
%       the triangle, if not, repeat a-c with the next triangle.  If it is
%       inside the triangle, that is the point on the tissue surface which
%       corresponds to the fiber.  Now we can move this inside the mesh
%       along the line segment connecting the 2 opposing fibers for any
%       distance (1/mus' for sources and just inside for detectors.


% Inputs:
% mesh is the mesh filename or workspace variable complete with mesh.source
% and mesh.meas information produced from mimics2nirfast_sd_3D_mouse
% w is width of gaussian source
% fn_oppfibers is a text file which specifies which fibers are directly
% opposite one another.  The fiber coordinates will be moved along the line
% connecting these positions.

if ischar(mesh)==1
    mesh = load_mesh(mesh);
end

% Make sure first column of opposing fibers list is 1:#sources, matching
% the source coord.
opp_fibers = load(fn_oppfibers);

%****************************************************
% Move opposing fibers toward one another and place on mesh surface:

[nsource,junk] = size(mesh.source.coord);
% error check
if length(opp_fibers)~=nsource
    disp('Error: length of opposing fibers list does not match length of source.coord')
    return
end

for i = 1 : nsource
    % If "guess" of s-d coordinate (the coordinates input into
    % this program) is inside mesh surface, plane_line_interesect will
    % return "Does not intersect".  So, extend line segment to
    % ensure it intersects with test plane:
    vec = (mesh.source.coord(opp_fibers(i,1),:)-mesh.source.coord(opp_fibers(i,2),:));
    uvec = vec./sqrt(sum(vec.^2));
    % Add 20 mm along segment on current s-d coordinate
    newcoord = 20*uvec+mesh.source.coord(opp_fibers(i,1),:);
    clear vec uvec
    
    % distance of each boundary node from source
    dist = distance(mesh.nodes,mesh.bndvtx,[mesh.source.coord(i,:) 0]);
    
    % sort nodes in order of nearest to farthest
    [snodes,ind] = sort(dist);
    
    % Start with nearest node and test surface elements of which node is a
    % vertex.  Call this node 'x'
    j=1; true = 0;
    % For each test node 'x'
    while true == 0
        % ind is list of node numbers in order of distance from
        % source/det point.  Find the elements which contain node 'x':
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
            % Find equation of line between two surface nodes
            A = [mesh.nodes(foo2(k,1),1),1; mesh.nodes(foo2(k,2),1),1];
            b = [mesh.nodes(foo2(k,1),2); mesh.nodes(foo2(k,2),2)];
            mb1 = A\b;
            
            % Find equation of line between fibers
            A = [mesh.source.coord(opp_fibers(i,1),1),1;mesh.source.coord(opp_fibers(i,2),1),1];
            b = [mesh.source.coord(opp_fibers(i,1),2);mesh.source.coord(opp_fibers(i,2),2)];
            mb2 = A\b;
            
            % Find intersection point of two lines:
            A = [mb1(1),-1; mb2(1),-1];
            b = [-mb1(2);-mb2(2)];
            xy = A\b;
            
            if (isnan(xy)) == 0
                true = 1;
                mesh.source.coord(opp_fibers(i,1),:) = xy;
            else
                true = 0;
                k = k+1;
            end
        end
        j=j+1;
    end
end
mesh.meas.coord = mesh.source.coord;

%***********************************************
% Now move sources in 1 scatter distance from surface and detectors just
% inside mesh (by 0.001)
if strcmp(mesh.type,'fluor')==1
    mus_eff = mesh.musx;
elseif strcmp(mesh.type,'stnd')==1
    mus_eff = mesh.mus;
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
[ind,int_func] = mytsearchn(mesh,...
    mesh.meas.coord);
mesh.meas.int_func = [ind int_func]


%*******************************************
% Finally, fix calculated positions so s-d pairs are not moved when loaded
% using load_mesh
mesh.meas.fixed = 1;
mesh.source.fixed = 1;



% Subfunction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,check]=plane_line_intersect(n,V0,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs:
%       n: normal vector of the Plane
%       V0: any point that belongs to the Plane
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
%Outputs:
%      I    is the point of interection
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1
%
% Example:
% Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
% The plane is represented by the normal vector n=[1 1 1]
% and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
% The segment is represented by the following two points
% P0=[-5 1 -1]
%P1=[1 2 3]
% [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);

%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.

I=[0 0 0];
u = P1-P0;
w = P0 - V0;
D = dot(n,u);
N = -dot(n,w);
check=0;
if abs(D) < 10^-7        % The segment is parallel to plane
    if N == 0           % The segment lies in plane
        check=2;
        return
    else
        check=0;       %no intersection
        return
    end
end

%compute the intersection parameter
sI = N / D;
I = P0+ sI.*u;

if (sI < 0 || sI > 1)
    check= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    check=1;
end


