%%to find the data points above & below a certain threshold of voltage to
%%remove noisy data
%%created subha 4/10/07
%%thresh must have 2 elements, low and high thresholds.
function [link] = remove_meas_mike(fn_data,thresh)
d = importdata([fn_data(1:end-3) 'dat'],'\t',61);
vol = d.data(:,2);
link = vol > thresh(1) & vol < thresh(2);
% for looking at raw data
%load('C:\Users\michael_a_mastanduno\Documents\LINK changes\full_link.mat')
%test = [full_link(:,1:2) repmat((1:15)',16,1) vol link]
end