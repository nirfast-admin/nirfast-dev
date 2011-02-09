function [paa] = drop_dist(paa, mesh, dist_thresh)
% drop_dist drops data points below a certain distance threshold
%   paa must be a 240xN matrix of data
%   mesh is for calculating distances
%   dist thresh is in mm. usually 15.

dist = sqrt(sum((mesh.source.coord(paa.link(:,1),:) - ...
    mesh.meas.coord(paa.link(:,2),:)).^2,2));

a = dist<dist_thresh;
paa.link(a,3:end) = 0;
end
    
