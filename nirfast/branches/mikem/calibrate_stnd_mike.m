function [data,mesh] = calibrate_stnd_mike(homog_data,...
                   anom_data,...
                   mesh_homog,...
                   mesh_anom,...
        		   frequency,...
        		   iteration)

% [data,mesh] = calibrate_stnd(homog_data,...
%                   anom_data,...
%                   mesh_homog,...
%                   mesh_anom,...
%        		   frequency,...
%        		   iteration)
%
% Calibrates standard data, using homogeneous data measured on phantom
% and the anomaly data
%
% homog_data is the homogeneous data (filename or variable)
% anom_data is the anomaly data (filename or variable)
% mesh_homog is the homogeneous mesh (filename or variable)
% mesh_anom is the anomaly mesh (filename or variable)
% frequency is the modulation frequency (MHz)
% iteration is the number of iterations for fitting
% data is the calibrated data
% mesh is the calibrated mesh with initial guesses



% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% If not a workspace variable, load meshes
if ischar(mesh_homog)== 1
  mesh_homog = load_mesh(mesh_homog);
end

if ischar(mesh_anom)== 1
  mesh_anom = load_mesh(mesh_anom);
end

mesh = mesh_anom;

% load anomaly data
paa_anom   = load_data(anom_data);
if ~isfield(paa_anom,'paa')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end
paa_anom = paa_anom.paa;

% set phase in radians
[j,k] = size(paa_anom);
% for i=1:2:k
% paa_anom(find(paa_anom(:,i+1)<0),i+1) = ...
%     paa_anom(find(paa_anom(:,i+1)<0),i+1) + (360);
% paa_anom(find(paa_anom(:,i+1)>(360)),i+1) = ...
%     paa_anom(find(paa_anom(:,i+1)>(360)),i+1) - (360);
% end

% load homogeneous data
paa_homog = load_data(homog_data);
paa_homog = paa_homog.paa;

% set phase in radians
[j,k] = size(paa_homog);
% for i=1:2:k
% paa_homog(find(paa_homog(:,i+1)<0),i+1) = ...
%     paa_homog(find(paa_homog(:,i+1)<0),i+1) + (360);
% paa_homog(find(paa_homog(:,i+1)>(360)),i+1) = ...
%     paa_homog(find(paa_homog(:,i+1)>(360)),i+1) - (360);
% end


    
% Calculate global mua and mus plus offsets for phantom data
if frequency == 0
    [mua_h,mus_h,lnI_h,data_h_fem] = fit_data_cw(mesh_homog,...
                                                      paa_homog(:,1),...
                                                      iteration);
    phase_h = zeros(size(lnI_h,1),1);
else
    [mua_h,mus_h,lnI_h,phase_h,data_h_fem] = fit_data_mike(mesh_homog,...
                                                  paa_homog,...
                                                  frequency,...
                                                  iteration);
end

% Calculate global mua and mus plus offsets for patient data
if frequency == 0
    [mua_a,mus_a,lnI_a,data_a_fem] = fit_data_cw(mesh_anom,...
                                                      paa_anom(:,1),...
                                                      iteration);
    phase_a = zeros(size(lnI_a,1),1);
else
    [mua_a,mus_a,lnI_a,phase_a,data_a_fem] = fit_data_mike(mesh_anom,...
                                                  paa_anom,...
                                                  frequency,...
                                                  iteration);
end

% calculate offsets between modeled homogeneous and measured
% homogenous and using these calibrate data
data_h_fem(:,1) = log(data_h_fem(:,1));
data_a_fem(:,1) = log(data_a_fem(:,1));
paa_anom(:,1) = log(paa_anom(:,1));
paa_homog(:,1) = log(paa_homog(:,1));
% paa_anomtmp = paa_anom;
% paa_homogtmp = paa_homog;
paa_cal = paa_anom - ((paa_homog - data_h_fem));


% %%%%%%%%%%%  Here are the lines in question...  %%%%%%%%%%%%%%
%paa_cal(:,1) = paa_cal(:,1) - (lnI_a-lnI_h);
%paa_cal(:,2) = paa_cal(:,2) - (phase_a-phase_h);

% This block applies the reflection intercept to all data
[data_a_fem_amp_line data_a_fem_phase_line] = plot_lnri_ref(data_a_fem,mesh_anom);
[paa_cal_amp_line paa_cal_phase_line] = plot_lnri_ref(paa_cal,mesh_anom);
paa_cal(:,1) = paa_cal(:,1) - (paa_cal_amp_line(2) - data_a_fem_amp_line(2));
paa_cal(:,2) = paa_cal(:,2) - (paa_cal_phase_line(2) - data_a_fem_phase_line(2));

% This block applies reflection and transmission intercepts to each data
% seperately.
% ref = importdata('link_reflectance.txt');
% ref = reshape(ref',240,1);
% ref = ref~=0;
% trans = ~ref;
% % reflection points
% [data_a_fem_amp_line_ref data_a_fem_phase_line_ref] = plot_lnri_ref(data_a_fem,mesh_anom);
% [paa_cal_amp_line_ref paa_cal_phase_line_ref] = plot_lnri_ref(paa_cal,mesh_anom);
% paa_cal(ref,1) = paa_cal(ref,1) - (paa_cal_amp_line_ref(2) - data_a_fem_amp_line_ref(2));
% paa_cal(ref,2) = paa_cal(ref,2) - (paa_cal_phase_line_ref(2) - data_a_fem_phase_line_ref(2));
% % transmission points
% [data_a_fem_amp_line_trans data_a_fem_phase_line_trans] = plot_lnri_trans(data_a_fem,mesh_anom);
% [paa_cal_amp_line_trans paa_cal_phase_line_trans] = plot_lnri_trans(paa_cal,mesh_anom);
% paa_cal(trans,1) = paa_cal(trans,1) - (paa_cal_amp_line_trans(2) - data_a_fem_amp_line_trans(2));
% paa_cal(trans,2) = paa_cal(trans,2) - (paa_cal_phase_line_trans(2) - data_a_fem_phase_line_trans(2));
% % %%%%%%%%%

paa_cal(:,1) = exp(paa_cal(:,1));
% calibrated data out
data.paa = paa_cal;
   
if 1 == 1
    % make plots of amplitude and homog data on top of one another 
    AH = paa_homog(:,1);
    PH = paa_homog(:,2);
    AA = paa_anom(:,1);
    PA = paa_anom(:,2);
    hnan = isnan(AH);
    anan = isnan(AA);
    AH(anan) = NaN;
    PH(anan) = NaN;
    AA(hnan) = NaN;
    PA(hnan) = NaN;
    anan = isnan(AH);
    AH = AH(~anan);
    PH = PH(~anan);
    AA = AA(~anan);
    PA = PA(~anan);
    figure;
    subplot(2,2,1)
    plot(AH,'b')
    hold on;
    plot(AA,'r')
    title('Amplitude')
    legend('Homog','Anom')
    subplot(2,2,2)
    plot(PH,'b')
    hold on;
    plot(PA,'r')
    title('Phase')
    legend('Homog','Anom')
    % plot calibrated data and simulated data on anomaly mesh
    Afem = data_a_fem(:,1);
    Pfem = data_a_fem(:,2);
    Acal = paa_cal(:,1);
    Pcal = paa_cal(:,2);
    femnan = isnan(Afem);
    calnan = isnan(Acal);
    Afem(calnan) = NaN;
    Pfem(calnan) = NaN;
    Acal(femnan) = NaN;
    Pcal(femnan) = NaN;
    calnan = isnan(Afem);
    Afem = Afem(~calnan);
    Pfem = Pfem(~calnan);
    Acal = Acal(~calnan);
    Pcal = Pcal(~calnan);
    subplot(2,2,3)
    plot(Afem,'b')
    hold on;
    plot(log(Acal),'r')
    title('Amplitude')
    legend('data a fem','paa cal')
    subplot(2,2,4)
    plot(Pfem,'b')
    hold on;
    plot(Pcal,'r')
    title('Phase')
    legend('data a fem','paa cal')
end
% set mesh values for global calculated patient values
mesh.mua(:) = mua_a;
mesh.mus(:) = mus_a;
mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
    