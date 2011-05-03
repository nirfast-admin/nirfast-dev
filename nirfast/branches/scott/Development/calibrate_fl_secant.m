function [data,mesh] = calibrate_fl_secant(fmesh, data_meas,...
    frequency, iteration, tolerance)

% [data,mesh] = calibrate_fl(fmesh, data_meas, ...
%    frequency, iteration, tolerance)
%
% Calibrates fluor data and generates initial guess. Will skip
% the calibrate part if only fluorescence data is available.
%
% fmesh is the forward mesh (variable or filename)
% data_meas is the uncalibrated data (variable or filename)
%       contains data_meas.amplitudex and data_meas.amplitudefl
% frequency is the modulation frequency (MHz)
% iteration is the number of iterations for fitting
% tolerance is the fitting tolerance to stop at
% data is the resulting calibrated data
% mesh is the resulting calibrated mesh with initial guess


% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

if ~exist('tolerance','var')
    tolerance = 10e-5;
end

% If not a workspace variable, load mesh
if ischar(fmesh)== 1
    fmesh = load_mesh(fmesh);
end
mesh = fmesh; clear fmesh

% load data
if ischar(data_meas)
    data_meas = load_data(data_meas);
end

if ~isfield(data_meas,'amplitudefl') && ~isfield(data_meas,'link')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end

% ***********************************************************
% Calibrate data to model:

% run forward model to get modeled excitation field
disp('Calculating excitation field...')
mesh.link = data_meas.link;
data_fwd = femdata(mesh, frequency);

% calculate calibrated data
if isfield(data_meas,'amplitudex') && isfield(data_fwd,'amplitudex')
    data.amplitudefl = data_meas.amplitudefl.*(data_fwd.amplitudex./data_meas.amplitudex);
else
    data.amplitudefl = data_meas.amplitudefl;
end
data.link = data_meas.link;
clear data_meas

ind = data.link(:,3)==0;
tempdata = data.amplitudefl;
tempdata(ind,:) = [];

% ***********************************************************
% Use the secant method to generate initial homogeneous guess
disp('Initializing Secant method points...')
lnI = log(tempdata);
muaf0 = 10^-7;
muaf1 = 5^-7;

mesh.muaf(:) = muaf0;
Err0 = femdata_err(mesh,frequency,lnI,ind);

mesh.muaf(:) = muaf1;
Err1 = femdata_err(mesh,frequency,lnI,ind);

for i = 1:iteration
    muaf2 = muaf1 - Err1*(muaf1-muaf0)/(Err1-Err0);
    muaf0 = muaf1;
    Err0 = Err1;
    muaf1 = muaf2;
    mesh.muaf(:) = muaf1;
    Err1 = femdata_err(mesh,frequency,lnI,ind);
    
    if i>1 && abs(Err1-Err0) < tolerance
        disp(['Stopping Criteria Reached at iteration ' num2str(i)]);
        disp('Global values calculated from Numerical fit');
        disp(['muaf = ' num2str(muaf1) ' mm-1 with error of ' num2str(Err1)]);
        disp('-------------------------------------------------');
        return
    end
    disp(['Iteration = ' num2str(i)]);
    disp('Global values calculated from Numerical fit');
    disp(['muaf = ' num2str(muaf1) ' mm-1 with error of ' num2str(Err1)]);
    disp('-------------------------------------------------');
end


function Err = femdata_err(mesh,frequency,lndata,ind)
[fem_data]=femdata(mesh,frequency);
fem_data.amplitudefl(ind,:) = [];
fem_lnI = log(fem_data.amplitudefl);
Err = sum((fem_lnI-lndata).^2);