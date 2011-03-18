function [loginfo] = load_log_file(data_fn)

%  Assumes log file is formatted with text at the top and a matrix of
%  spectrometer/camera parameters.  The list 16 columns of this matrix are
%  the exposure times used where each column corresponds to a source
%  position and each row a detector. An exposure time = zero indicates that
%  source-detector pair was not used.

% s. c. davis 2007

logtemp = importdata([data_fn,'.log']);
loginfo.exp_times = logtemp.data(:,end-15:end);
loginfo.filters = logtemp.data(:,end-31:end-16);
loginfo.cameratemp = logtemp.data(:,3);
loginfo.centerwv = logtemp.data(:,7);
loginfo.sources = logtemp.data(:,1);
loginfo.meas = logtemp.data(:,2);
loginfo.grating = logtemp.data(:,6);
for i = 1:numel(loginfo.grating)
    if loginfo.grating(i) == 1
        loginfo.grating(i) = 1200;
    elseif loginfo.grating(i) == 2
        loginfo.grating(i) = 300;
    end
end

% create a nirfast-style link matrix based on sources and detectors
% selected in the acquisition program

temp_src = loginfo.sources;
temp_meas = loginfo.meas;

in0 = [];
for i = 1:numel(temp_src)
    if temp_src(i) == 0 & temp_meas(i) == 0
        in0 = [in0 i];
    end
end
temp_src(in0) = [];  temp_meas(in0) = [];
loginfo.sources(in0) = NaN; loginfo.meas(in0) = NaN; 


cs = 1; cd = 1;
ns = []; nd = [];
for i = 1:numel(temp_src)
    if temp_src(i) == 1 & temp_meas(i) == 1
        ns = [ns,cs]; nd = [nd,cd];
    elseif temp_src(i) == 1 & temp_meas(i) == 0
        ns = [ns,cs]; nd = [nd,0];
    elseif temp_src(i) == 0 & temp_meas(i) == 1
        ns = [ns,0]; nd = [nd,cd];
    end
    cs = cs+1; cd = cd+1;
end

% biuld link file based on sd_pairs
sd_ind = 1:length(ns);
cn = 1;
for i = 1:length(ns)
    if ns(i) ~= 0
        loginfo.link(cn,:) = nd([(i+1):end, 1:(i-1)]);
        cn = cn+1;
    end
end
%loginfo.link;
