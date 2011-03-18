function [m0 m1] = plot_lnri_cw(paa,mesh)
% plot_lnri Plots a data set as ln(rI) vs. sd distances the way calibrate
% does it. Alternative way to view data from analyze_data and view_data
%   paa is a 2 column data file to be plotted
%   mesh is the mesh which the SD distances will be drawn from for plotting
%   wv is the wavelength of the paa for the plot title
% fix phase wrap as it is done in calibrate
data = paa.paa;
[j,k] = size(data);
for i=1:k
clear m0 m1
data = paa.paa(:,i);
% get an index from link file of data to actually use
linki = logical(paa.link(:,i+2));
% calculate the source / detector distance for each combination.
dist = sqrt(sum((mesh.source.coord(paa.link(:,1),:) - ...
    mesh.meas.coord(paa.link(:,2),:)).^2,2));

% Set lnrI, lnr and phase!
lnrI = log(data(:,1).*dist);
lnI = log(data(:,1));

figure;
plot(dist(linki),lnrI(linki),'.')
ylabel('lnrI');
xlabel('Source / Detector distance');
drawnow
pause(0.001)

% Calculate the coeff of a polynomial fit of distance vs. Phase or lnrI
% then add fit lines to graph
m1 = polyfit(dist(linki),lnrI(linki),1);
x = min(dist(linki)):(max(dist(linki))-min(dist(linki)))/10:max(dist(linki));
hold on
plot(x,m1(1).*x+m1(2))
drawnow
axis square
pause(0.001)
end
end

