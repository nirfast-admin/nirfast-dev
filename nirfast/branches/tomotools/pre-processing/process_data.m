function [ data, mesh , dist] = process_data(mesh, datafile, cutoff)
%Kelly E. Michaelsen - 09/30/11

mesh=load_mesh(mesh);
data= importdata(datafile);

sources=11*7;
threshold=[0.1,9.5];
starttime=800;
wv=mesh.wv;
edge_limit=20;
%wv=mesh.wv(1:4);

%The following functions convert the data from the labview file to a
%nirfast compatible data file

data=sd_average(data, starttime);
data=drop_data(data, threshold);
data=wv_data(data, wv);
data=detector_calib(data);
%data=source_calib(data, sources); %Only good for the 20 source geometry
data=data_for_nirfast(data,sources, wv);

%The following functions perform data preprocessing by selectively removing
%data points based on source and detector positions as well as distances

mesh.link=data.link;
%Move SD from original position to optimize data
mesh.source.coord(:,2)=mesh.source.coord(:,2)+1;
mesh.source.coord(:,1)=mesh.source.coord(:,1)-1;
[data,dist]=data_rm_far(data, mesh, cutoff);
[data]=data_rm_edges(data,mesh);
[data]=data_rm_source_edges(data,mesh,edge_limit);
%[data]=data_rm_bottle_detectors(data, mesh);
%[data]=data_rm_62(data);

%Plot the results
plot_data(data, dist);
end


function [ data ] = sd_average( data, starttime )
%This function takes the data as well as a starting measurement
%All measurements before this starting measurement are discarded due to
%laser instability after being switched on
%An average measurement of the leftover measurements is taken for each SD
%pair
data=data(:,starttime:end);
data=mean(data,2);
end

function [ data ] = drop_data( data, threshold )
%This function drops all the data points below a certain threshold (at the
%noise floor) and above another threshold (saturation)
data(data<threshold(1))=NaN;
data(data>threshold(2))=NaN;
end

function [ wv_data ] = wv_data( data, wv)
%This function creates a column for each wavelength of data
count=1;
data_per_wv=length(data)/length(wv);
for i=1:length(wv)
    wv_data(:,i)=data(count:count+data_per_wv-1);
    count=count+data_per_wv;
end
end

function [ organized_data ] = detector_calib( data)
%This function sorts the data so that it is in the correct detector order
%(detector 1 corresponds to the back left corner, detector two is to the
%right of detector 1 across the row, then restarting at the left corner of
%the second row etc...)
%Also multiplies the data by a sensitivity factor determined from detector
%calibration for each detector
detect_LUT=load ('LUT_detector.txt');
organized_data=[];
[datapoints, numb_wv]=size(data);
for i=1:numb_wv
    organized_data_wv=[];
    for count=1:80:datapoints
        data_single_source=data(count:count+79,i);
        org_single_source=data_single_source(detect_LUT(:,1));
        org_single_source=org_single_source./detect_LUT(:,2);
        organized_data_wv=[organized_data_wv;org_single_source];
    end
    organized_data(:,i)=organized_data_wv;
end
end

function [ calib_data ] = source_calib(data, numb_sources)
%This function alters the data by a sensitivity factor determined via
%calibration for each wavelength and source position only good for the
%original 20 source setup
calib_data=[];
source_LUT=load ('LUT_source.txt');
[datapoints, numb_wv]=size(data);
for i=1:numb_wv
    calib_data_wv=[];
    source_numb=1;
    for count=1:75:datapoints
        data_single_source=data(count:count+74,i);
        calib_data_single_source=data_single_source./source_LUT(source_numb,i);
        calib_data_wv=[calib_data_wv; calib_data_single_source];
        source_numb=source_numb+1;
    end
    calib_data(:,i)=calib_data_wv;
end
end

function [ nirfast_data ] = data_for_nirfast(data,numb_sources, wv)
%This function creates a nirfast compatible spectral datafile
calib_data=[];
source_LUT=load ('LUT_source.txt');
[datapoints, numb_wv]=size(data);
nirfast_data.link=ones(datapoints, numb_wv+2);
for ind=1:datapoints
    nirfast_data.link(ind,2)=mod(ind-1,75)+1;
    nirfast_data.link(ind,1)=floor(abs(ind-1)/75)+1;
end
nirfast_data.paa=zeros(datapoints,2*numb_wv);
nirfast_data.paa(:,1:2:numb_wv*2)=data;
nirfast_data.wv=wv;
end

function [ data, dist_full ] = data_rm_far(data,mesh, cutoff)
%This function removes data points with source detector distances that are 
%further than a cutoff distance
mesh.link=data.link;
for i = 1:length(mesh.link)
    snum = mesh.link(i,1);
    mnum = mesh.link(i,2);
    snum = mesh.source.num == snum;
    mnum = mesh.meas.num == mnum;
    if sum(snum)==0 || sum(mnum)==0
        dist_full(i,1)=0;
        mesh.link(i,3)=0;
    else
        dist_full(i,1) = sqrt(sum((mesh.source.coord(snum,:) - ...
            mesh.meas.coord(mnum,:)).^2,2));
    end
end
far=find(dist_full>cutoff);
data.paa(far,:)=NaN;
dist_full(far)=NaN;
end
function [ data ] = data_rm_edges(data,mesh)
%This function removes data points from the exterior most detectors
xedge=[max(mesh.meas.coord(:,1)),min(mesh.meas.coord(:,1))];
yedge=[max(mesh.meas.coord(:,2)),min(mesh.meas.coord(:,2))];
numb_detect=length(mesh.meas.coord);
edge_detector=find(mesh.meas.coord(:,1)==xedge(1) | mesh.meas.coord(:,1)==xedge(2)...
    |mesh.meas.coord(:,2)==yedge(1) | mesh.meas.coord(:,2)==yedge(2));
for i=1:length(edge_detector)
    data.paa(edge_detector(i):numb_detect:end,:)=NaN;
end
end

function [data]=data_rm_source_edges(data,mesh,lim)
%This function removes data points for any sources less than a cutoff 
%distance from the edge of the mesh
xedge=[max(mesh.nodes(:,1)),min(mesh.nodes(:,1))];
yedge=[max(mesh.nodes(:,2)),min(mesh.nodes(:,2))];
numb_detect=length(mesh.meas.coord);
edge_source=[find(abs(mesh.source.coord(:,1)-xedge(1))<lim); find(abs(mesh.source.coord(:,1)-xedge(2))<lim);...
    find(abs(mesh.source.coord(:,2)-yedge(1))<lim); find(abs(mesh.source.coord(:,2)-yedge(2))<lim)];
for i=1:length(edge_source)
    data.paa(numb_detect*(edge_source(i)-1)+1:numb_detect*(edge_source(i)),:)=NaN;
end
end

function plot_data(data, dist)
    figure
    for i=1:length(data.wv)
    subplot(2, length(data.wv)/2,i)
    wv_data=data.paa(:,(i-1)*2+1);
    plot(dist,log(wv_data),'.','MarkerSize', 15 )
    set(gcf, 'color', 'white');
    xlabel('Source-Detector Distance in mm', 'FontSize', 18)
    ylabel('Ln Amplitude of Calibrated Detector Response', 'FontSize', 18)
    title(['Ln Amplitude at ', num2str(data.wv(i)),' nm vs SD Distance'], 'FontSize', 20)
    end
end

% 
% function [ data ] = data_rm_bottle_detectors(data,mesh, numb_sources)
% %This function removes detectors out of range of measurements of bottle...
% xedge=[max(mesh.meas.coord(:,1)),min(mesh.meas.coord(:,1))];
% yedge=[max(mesh.meas.coord(:,2)),min(mesh.meas.coord(:,2))];
% numb_detect=length(mesh.meas.coord);
% edge_detector=find(mesh.meas.coord(:,1)==xedge(1) | mesh.meas.coord(:,1)==xedge(2)...
%     |mesh.meas.coord(:,2)==yedge(1) | mesh.meas.coord(:,2)==yedge(2));
% for i=1:length(edge_detector)
%     data.paa(edge_detector(i):numb_detect:end,:)=NaN;
% end
% end
