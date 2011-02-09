function [adata hdata amesh hmesh] = dataprep(a_fn,h_fn,amesh,hmesh,arep,hrep,aplane,hplane,fibers)
% [hdata adata] = dataprep(h_fn,a_fn,hmesh,amesh)

% Drops noisy data based on PMT S/N
% Drops near source data and bad fibers
% Data files to use are adata and hdata
close all
load('C:\Users\michael_a_mastanduno\Documents\LINK changes\full_link.mat')
link = full_link;
% if plotflag is set to 1, it will plot data lnri vs. sd distance
plotflag = 1;
wavelengths = [661 735 785 808 826 849];

drop_thresh = [.05 4];
dist_thresh = 16;

% HOMOG load and drop data, then convert to spectral
disp(['Dropping HOMOG data at ',num2str(wavelengths)])
i = 3;
j = 1;
for wv = wavelengths
    hpaa = [num2str(h_fn),num2str(wv),'nm_rep',num2str(hrep),'_plane',num2str(hplane),'.paa'];
    if fopen(hpaa) == -1
        hpaa = [num2str(h_fn),num2str(wv),'nm_rep',num2str(hrep),'.paa'];
    end
    tmp = load_data(hpaa);
    hdata.paa(:,j:j+1) = tmp.paa;
    link(:,i) = remove_meas_mike(hpaa,drop_thresh);
    i = i+1;
    j = j+2;
end
hdata.wv = wavelengths;
hdata.link = link;
% drop bad fibers
[hdata.link] = drop_fibers(hdata.link,fibers);
% drop based on distance
hdata = drop_dist(hdata,hmesh,dist_thresh);
hmesh.link = hdata.link;
% plot each wv for inspection.
if plotflag==1
    plot_lnri(hdata,hmesh);
end
i = 3;
for wv = wavelengths
    disp(['Using ',num2str(sum(link(:,i))),' data points at ', num2str(wv),'nm'])
    i = i + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANOMALY load and drop data, then convert to spectral
link = full_link;
disp(' ')
disp(['Dropping ANOM data at ',num2str(wavelengths)])
i = 3;
j = 1;
for wv = wavelengths
    apaa = [num2str(a_fn),num2str(wv),'nm_rep',num2str(arep),'_plane',num2str(aplane),'.paa'];
    if fopen(apaa) == -1
        apaa = [num2str(a_fn),num2str(wv),'nm_rep',num2str(arep),'.paa'];
    end
    tmp = load_data(apaa);
    adata.paa(:,j:j+1) = tmp.paa;
    link(:,i) = remove_meas_mike(apaa,drop_thresh);
    i = i+1;
    j = j+2;
end
adata.wv = wavelengths;
adata.link = link;
% drop bad fibers
[adata.link] = drop_fibers(adata.link,fibers);
% drop based on distance
adata = drop_dist(adata,amesh,dist_thresh);
amesh.link = adata.link;
% plot each wv for inspection.
if plotflag==1
    plot_lnri(adata,amesh);
end
i = 3;
for wv = wavelengths
    disp(['Using ',num2str(sum(link(:,i))),' data points at ', num2str(wv),'nm'])
    i = i + 1;
end

clear a_fn apaa arep data drop_thresh h_fn hpaa hrep i wavelengths m0 m1
clear dist_thresh plotflag aplane hplane
end
