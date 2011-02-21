function [m0 m1] = plot_lnri(paa,mesh)
% plot_lnri Plots a data set as ln(rI) vs. sd distances the way calibrate
% does it. Alternative way to view data from analyze_data and view_data
%   paa is a 2 column data file to be plotted
%   mesh is the mesh which the SD distances will be drawn from for plotting
%   wv is the wavelength of the paa for the plot title
% fix phase wrap as it is done in calibrate
data = paa.paa;
[j,k] = size(data);
j = 1;
for i=1:k/2
clear m0 m1
data = paa.paa(:,j:j+1);
% get an index from link file of data to actually use
linki = logical(mesh.link(:,i+2));
% calculate the source / detector distance for each combination.
dist = sqrt(sum((mesh.source.coord(mesh.link(:,1),:) - ...
    mesh.meas.coord(mesh.link(:,2),:)).^2,2));

% Set lnrI, lnr and phase!
lnrI = log(data(:,1).*dist);
lnI = log(data(:,1));
phase = data(:,2);

figure;
subplot(1,2,1);
plot(dist(linki),lnrI(linki),'.')
ylabel('lnrI');
xlabel('Source / Detector distance');
subplot(1,2,2);
plot(dist(linki),phase(linki),'.')
ylabel('Phase');
xlabel('Source / Detector distance');
drawnow
pause(0.001)

% Calculate the coeff of a polynomial fit of distance vs. Phase or lnrI
% then add fit lines to graph
m0 = polyfit(dist(linki),phase(linki),1);
m1 = polyfit(dist(linki),lnrI(linki),1);
x = min(dist(linki)):(max(dist(linki))-min(dist(linki)))/10:max(dist(linki));
subplot(1,2,2)
hold on
plot(x,m0(1).*x+m0(2))
drawnow
pause(0.001)
axis square
subplot(1,2,1)
hold on
plot(x,m1(1).*x+m1(2))
drawnow
axis square
pause(0.001)
j=j+2;
end
end
