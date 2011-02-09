function [lnrI_line phase_line] = plot_lnri_trans(paa,mesh)
% [amp_slope amp_off phs_slope phs_off] = plot_lnri_ref(paa,mesh)
%
% plot_lnri Plots a data set as ln(rI) vs. sd distances the way calibrate
% does it. Alternative way to view data from analyze_data and view_data
%   paa is a 2 column data file to be plotted
%   mesh is the mesh which the SD distances will be drawn from for plotting
%   wv is the wavelength of the paa for the plot title
% fix phase wrap as it is done in calibrate

[j,k] = size(paa);
link_orig = mesh.link;
clear m0 m1
data = paa;
if data(~isnan(data(:,1)),1)<0
    data(:,1) = exp(data(:,1));
end
link_t = importdata('link_transmission.txt');
mesh.link = link_t;
a = link_orig==0;
mesh.link(a) = 0;
% calculate the source / detector distance for each measurement
k = 1;
datanum = 0;
[ns,junk]=size(mesh.source.coord);
for i = 1 : ns
    for j = 1 : length(mesh.link(i,:))
        datanum = datanum + 1;
        if isnan(data(datanum,1)) || isnan(data(datanum,2))
            mesh.link(i,j) = 0;
        end
        if mesh.link(i,j) ~= 0
            jj = mesh.link(i,j);
            dist(k,1) = sqrt(sum((mesh.source.coord(i,:) - ...
                mesh.meas.coord(jj,:)).^2));
            k = k+1;
        else
            data(datanum,1:2) = NaN;
        end
    end
end

% Set lnrI, lnr and phase!
[j,k] = size(data(:,1));
[j2,k2] = size(dist);

% deal with NaNs
dist_orig = dist;
ind = unique([find(isnan(data(:,1))==1); find(isnan(data(:,2))==1)]);
ind = setdiff(1:size(data,1),ind);
data = data(ind,:);



% multiply for plotting
lnrI = log(data(:,k).*dist);
lnI = log(data(:,k));
phase = data(:,k+1);

% figure;
% subplot(1,2,1);
% plot(dist,lnrI,'b.')
% ylabel('lnrI');
% xlabel('Source / Detector distance');
% subplot(1,2,2);
% plot(dist,phase,'b.')
% ylabel('Phase');
% xlabel('Source / Detector distance');
% drawnow
% pause(0.001)

% Calculate the coeff of a polynomial fit of distance vs. Phase or lnrI
% then add fit lines to graph
m0 = polyfit(dist,phase,1);
m1 = polyfit(dist,lnrI,1);
% x = min(dist):(max(dist)-min(dist))/10:max(dist);
% subplot(1,2,2)
% hold on
% plot(x,m0(1).*x+m0(2))
% drawnow
% pause(0.001)
% axis square
% subplot(1,2,1)
% hold on
% plot(x,m1(1).*x+m1(2))
% drawnow
% axis square
% pause(0.001)

lnrI_line = m1;
phase_line = m0;
end

