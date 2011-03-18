function [data] = remove_cursor_data_func(ci,data,mesh,wv_num)
% 1. Export cursor data from a figure named by wavelength. 
% 2. Run this function to set data.link to 0 in that place.
% Assumes 6 wavelength data set.


% get an index from link file of data to actually use
linki = logical(data.link(:,wv_num+2));
% calculate the source / detector distance for each combination.
dist = sqrt(sum((mesh.source.coord(data.link(:,1),:) - ...
    mesh.meas.coord(data.link(:,2),:)).^2,2));
lnrI = log(data.paa(:,wv_num*2-1).*dist);
for i = 1:length(ci)
    % amplitude points
    a1 = and(dist>=ci(1,i).Position(1)-0.001,...
        dist<=ci(1,i).Position(1)+0.001);
    a2 = and(lnrI>=ci(1,i).Position(2)-0.001,...
        lnrI<=ci(1,i).Position(2)+0.001);
    a = and(a1,a2);
    % phase points
    b1 = and(dist>=ci(1,i).Position(1)-0.001,...
        dist<=ci(1,i).Position(1)+0.001);
    b2 = and(data.paa(:,wv_num*2)>=ci(1,i).Position(2)-0.001,...
        data.paa(:,wv_num*2)<=ci(1,i).Position(2)+0.001);
    b = and(b1,b2);
    c = or(a,b);
    nix(c) = 1;
end
nix = logical(nix);
data.link(nix,wv_num+2) = 0;
end

